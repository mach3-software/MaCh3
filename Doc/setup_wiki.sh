#Modules in case one neeeds them
#dnf install -y texlive-scheme-basic texlive-latex
#dnf install -y mathjax
#dnf install -y perl
#dnf install -y doxygen
#Download MaCh3 wiki
git clone https://github.com/mach3-software/MaCh3.wiki.git Doc/User-Guide
#Download fancy animation
BASE_URL="https://raw.githubusercontent.com/KSkwarczynski/Python/main/MaCh3/animation/"
curl -L "$BASE_URL/MCMC_Example.mp4" -o Doc/MaCh3Web/MCMC_Example.mp4
curl -L "$BASE_URL/posterior_predictive_Example.mp4" -o Doc/MaCh3Web/posterior_predictive_Example.mp4
curl -L "$BASE_URL/umbrella_sampling_delta_cp.mp4" -o Doc/MaCh3Web/umbrella_sampling_delta_cp.mp4

# Count publications and theses from the wiki
RESULTS_FILE="Doc/User-Guide/14.-MaCh3-in-the-Field.md"

PUBLICATIONS=$(awk '
    /^# Publications using MaCh3/ { found=1; next }
    /^# Theses using MaCh3/       { found=0 }
    found && /^\* /               { count++ }
    END { print count+0 }
' "$RESULTS_FILE")

THESES=$(awk '
    /^# Theses using MaCh3/ { found=1; next }
    found && /^\* /          { count++ }
    END { print count+0 }
' "$RESULTS_FILE")

echo "Found $PUBLICATIONS publications and $THESES theses"

# Update the numbers in the Doxygen source
sed -i -E "s/[0-9]+(\+)? papers using MaCh3/${PUBLICATIONS}+ papers using MaCh3/" Doc/MaCh3Web/mainpage.md
sed -i -E "s/[0-9]+(\+)? theses and dissertations/${THESES}+ theses and dissertations/" Doc/MaCh3Web/mainpage.md
