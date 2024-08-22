#ifndef STYLE_MANAGER_H
#define STYLE_MANAGER_H 1

#include "toml/toml_helper.h"

// ROOT Things
#include "TH1.h"
#include "TStyle.h"
#include "TColor.h"

namespace MaCh3Plotting{
    class StyleManager{
        public:

        StyleManager(std::string tomlName);

        std::string prettifyParamName(std::string origName) const;

        std::string prettifySampleName(std::string fullName) const;

        // style setting options
        void setPalette(int rootPlotStyle) const;
        void setPalette(std::string configStyleName) const;
        void setTH1Style(TH1* hist, std::string styleName) const;

        std::string prettifySampleName(std::string fullName);


        private:
        toml::value style_toml;
    };

}




#endif