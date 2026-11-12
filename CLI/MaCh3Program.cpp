/// @file MaCh3Program.cpp
/// @brief Implementation of the MaCh3Program class

#include <dlfcn.h>
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#include "CLI/MaCh3Program.hpp"

using std::vector, std::string, std::map;
namespace fs = std::filesystem;

namespace {
    // Fix executable paths: resolve natively on Linux/macOS, retaining a PATH fallback.
    std::string executable_path() {
        std::error_code error;
#ifdef __APPLE__
        std::uint32_t size = 1024;
        std::vector<char> buffer(size);
        if (_NSGetExecutablePath(buffer.data(), &size) != 0) {
            buffer.resize(size);
            if (_NSGetExecutablePath(buffer.data(), &size) != 0) return "mach3";
        }
        const auto path = fs::canonical(buffer.data(), error);
#else
        const auto path = fs::read_symlink("/proc/self/exe", error);
#endif
        return error ? "mach3" : path.string();
    }

    // Fix executable paths: prevent spaces, quotes or metacharacters changing shell parsing.
    std::string shell_quote(const std::string& value) {
        std::string result = "'";
        for (const char character : value) {
            if (character == '\'') result += "'\\''";
            else result += character;
        }
        return result + "'";
    }
}

namespace M3{
    // Fix completion installation: check open/write/close and report the failing operation.
    bool MaCh3Program::write_file(const fs::path& path, std::string_view content) const{
        std::string_view operation = "creating parent directories";
        try {
            errno = 0;
            fs::create_directories(path.parent_path());
            std::ofstream out;
            out.exceptions(std::ios::failbit | std::ios::badbit);
            operation = "opening";
            out.open(path);
            operation = "writing";
            out << content;
            operation = "closing";
            out.close();
            return true;
        } catch (const std::exception& err) {
            const int error = errno;
            std::cerr << "Cannot write completion script '" << path.string() << "' while " << operation
                      << ": " << err.what();
            if (error != 0) std::cerr << ": " << std::strerror(error);
            std::cerr << "\n";
            return false;
        }
    }

    std::string MaCh3Program::detect_shell() const{
        const char* shell = std::getenv("SHELL");
        if (!shell) return "";

        // Fix completion installation: recognise exact shell names, not substrings such as "notbash".
        const auto name = fs::path(shell).filename().string();
        if (name == "bash" || name == "zsh") return name;

        return "";
    }

    // Fix completion installation: environment and file failures must reach the caller as status 1.
    int MaCh3Program::install_completions() const{
        const char* home = std::getenv("HOME");
        if (!home || home[0] == '\0') {
            std::cerr << "Cannot determine HOME directory\n";
            return 1;
        }

        std::string shell = detect_shell();
        if (shell.empty()) {
            std::cerr << "Could not detect your shell (check it's exported). Supported: bash, zsh\n";
            return 1;
        }

        // Resolve the absolute path of this binary so the completion script
        // works even when mach3 is not on PATH.
        const std::string binary_path = shell_quote(executable_path());

        // Replace the generic "mach3" invocation in the template with the
        // absolute binary path so the script is self-contained.
        auto embed_path = [&](std::string_view tmpl) -> std::string {
            std::string s(tmpl);
            const std::string placeholder = "mach3 --complete";
            const std::string replacement = binary_path + " --complete";
            for (std::size_t pos = s.find(placeholder); pos != std::string::npos;
                 pos = s.find(placeholder, pos + replacement.size()))
                s.replace(pos, placeholder.size(), replacement);
            return s;
        };

        fs::path home_path(home);

        if (shell == "bash") {
            fs::path path = home_path / ".local/share/bash-completion/completions/mach3";
            if (!write_file(path, embed_path(MaCh3Program::BASH_COMPLETION))) return 1;
            std::cout << "Installed bash completions → " << path << "\n";
            return 0;
        }

        if (shell == "zsh") {
            fs::path path = home_path / ".local/share/zsh/site-functions/_mach3";
            if (!write_file(path, embed_path(MaCh3Program::ZSH_COMPLETION))) return 1;
            std::cout << "Installed zsh completions → " << path << "\n";
        }
        return 0;
    }

    void MaCh3Program::completions(const std::string& prefix, const std::vector<std::string>& context) const {
        for (const auto& name : get_completions(prefix, context))
            std::cout << name << "\n";
    }

    void MaCh3Program::parse_args(int argc, const char *const argv[]) {
        std::vector<std::string> args(argv, argv + argc);
        this->insert_default_subcommands(args);
        MaCh3ArgumentParser::parse_args(args);
    }

    // Fix registration safeguards: reject invalid parsers/collisions before mutating argparse's maps.
    void MaCh3Program::validate_parser(const MaCh3ArgumentParser* parser) const{
        if (!parser || parser->name().empty()) {
            throw std::invalid_argument("Module/plugin must provide a non-null parser with a non-empty command name");
        }
        if (m_subparser_map.count(parser->name()) != 0) {
            throw std::invalid_argument("Duplicate command name: " + parser->name());
        }
    }

    void MaCh3Program::add_core_module(IModule& module){
        // Fix registration safeguards: another get_parser() call could replace a registered parser.
        if (std::any_of(m_module_map.begin(), m_module_map.end(),
                        [&module](const auto& entry) { return entry.second == &module; })) {
            throw std::invalid_argument("Core module is already registered");
        }
        MaCh3ArgumentParser* parser = module.get_parser();
        validate_parser(parser);
        this->add_subparser(*parser);
        m_module_map[parser] = &module;
    }

    /// @brief Load dynamic plugins from paths specified in MACH3_PLUGINS environment variable
    ///
    /// Reads the MACH3_PLUGINS environment variable (colon-separated paths), loads each .so file,
    /// and registers the plugins with the program. Each plugin must export create_plugin and
    /// destroy_plugin functions.
    void MaCh3Program::load_dynamic_plugins(){
        const char* env = std::getenv("MACH3_PLUGINS");
        if (!env || env[0] == '\0') return;

        std::stringstream ss(env);
        std::string path;

        while (std::getline(ss, path, ':')) {
            if (path.empty()) {
                std::cerr << "Ignoring empty plugin path in MACH3_PLUGINS\n";
                continue;
            }
            std::vector<std::string> files;
            // Fix discovery errors: one inaccessible path must not block other optional plugins.
            try {
                files = this->expand_plugin_path(path);
            }
            catch (const fs::filesystem_error& err) {
                std::cerr << "Cannot discover plugins in '" << path << "': " << err.what() << "\n";
                continue;
            }
            for (const auto& sofile : files) {
                // Fix registration safeguards: canonical paths also deduplicate symlink aliases.
                if (!m_attempted_plugin_paths.insert(sofile).second) {
                    std::cerr << "Skipping duplicate plugin library: " << sofile << "\n";
                    continue;
                }
                try{
                    std::unique_ptr<DynamicPlugin> dlplugin = std::make_unique<DynamicPlugin>(sofile);
                    MaCh3ArgumentParser* parser = dlplugin->get_parser();
                    validate_parser(parser);
                    this->add_subparser(*parser);
                    m_dynamic_plugin_map[parser] = std::move(dlplugin);
                }
                catch (const std::exception& e) {
                    std::cerr << "Failed to load plugin from " << sofile << ": " << e.what() << "\n";
                    continue;
                }
                catch (...) {
                    std::cerr << "Failed to load plugin from " << sofile << ": unknown error\n";
                    continue;
                }
            }
        }
    }

    /// @brief Execute the selected subcommand
    ///
    /// Determines which subcommand was invoked and runs the corresponding module or plugin.
    ///
    /// @return Exit code from the executed module (0 on success)
    int MaCh3Program::Run(){
        // Dispatch to the top-level module; modules handle their own nested subcommands
        for (const auto& [parser, module] : m_module_map) {
            if (this->is_subcommand_used(*parser)) return module->Run();
        }
        for (const auto& [parser, dplugin] : m_dynamic_plugin_map) {
            if (this->is_subcommand_used(*parser)) return dplugin->Run();
        }
        return 0;
    }


    /// @brief Expand a plugin path to a list of shared library files
    ///
    /// If the path is a file, returns it directly. If it's a directory,
    /// returns supported shared libraries found in that directory.
    ///
    /// @param path Path to a plugin file or directory
    /// @return Vector of absolute paths to .so files
    std::vector<std::string> MaCh3Program::expand_plugin_path(const std::string& path) const {
        std::vector<std::string> result;

        fs::path p(path);

        if (!fs::exists(p)) {
            std::cerr << "Plugin path does not exist: " << path << "\n";
            return result;
        }

        if (fs::is_regular_file(p)) {
            // Fix registration safeguards: normalise relative paths and symlinks before deduplication.
            result.push_back(fs::canonical(p).string());
        }
        else if (fs::is_directory(p)) {
            // Fix platform coverage: also recognise native macOS shared-library suffixes.
            for (auto& entry : fs::directory_iterator(p)) {
                const auto extension = entry.path().extension();
                bool is_good_plugin = extension == ".so";
#ifdef __APPLE__
                is_good_plugin = is_good_plugin || extension == ".dylib";
#endif
                try{
                    is_good_plugin = is_good_plugin && entry.is_regular_file();
                }
                catch (const fs::filesystem_error& err) {
                    // Fix discovery errors: a broken candidate must not hide valid neighbours.
                    std::cerr << "Cannot inspect plugin '" << entry.path().string() << "': " << err.what() << "\n";
                    continue;
                }
                if (is_good_plugin) {
                    try{
                        result.push_back(fs::canonical(entry.path()).string());
                    }
                    catch (const fs::filesystem_error& err) {
                        std::cerr << "Cannot canonicalize plugin '" << entry.path().string() << "': " << err.what() << "\n";
                        continue;
                    }
                }
            }
        }
        else {
            std::cerr << "Invalid plugin path (not file or directory): " << path << "\n";
        }

        // Fix registration safeguards: make candidate precedence independent of directory iteration order.
        std::sort(result.begin(), result.end());
        return result;
    }
}
