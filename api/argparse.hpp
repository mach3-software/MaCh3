#pragma once
#include <argparse/argparse.hpp>

using argparse::ArgumentParser;

namespace mach3 {

    class MaCh3ArgumentParser: public ArgumentParser{
        public:
            using ArgumentParser::ArgumentParser;
            virtual ~MaCh3ArgumentParser() = default;

            const std::string name() const{
                return this->m_program_name;
            }
            const std::list<std::reference_wrapper<ArgumentParser>>& subparsers() const{
                return this->m_subparsers;
            }

            const MaCh3ArgumentParser& get_subcommand_used() const{
                for (const std::reference_wrapper<ArgumentParser>& subparser : this->m_subparsers) {
                    if (this->is_subcommand_used(subparser.get())) {
                        return static_cast<MaCh3ArgumentParser*>(&(subparser.get()))->get_subcommand_used();
                    }
                }
                return (*this);
            }
    };
};
