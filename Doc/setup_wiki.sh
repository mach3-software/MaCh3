#Modules in case one needs them
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

# ---- Most recent publication / PhD (first bullet of each list = newest) ----
MAINPAGE="Doc/MaCh3Web/mainpage.md"

# HTML-escape
prep() {
    printf '%s' "$1" \
      | sed -e 's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g' \
      | sed -e 's/[&|\\]/\\&/g'
}

# First publication bullet
LATEST_PUB=$(awk '
    /^# Publications using MaCh3/ { found=1; next }
    /^# Theses using MaCh3/       { found=0 }
    found && /^\* /               { sub(/^\* /, ""); print; exit }
' "$RESULTS_FILE")

# Drop trailing link etc
LATEST_PUB=$(printf '%s' "$LATEST_PUB" | sed -E '
    s/ \[https?:\/\/.*$//
    s/ (\*|arXiv:).*$//
    s/\*//g
    s/\.?[[:space:]]*$/./
')

# First thesis bullet: "* Name, [Title](url), (**YEAR**)"  ->  "Name — Title"
LATEST_PHD=$(awk '
    /^# Theses using MaCh3/ { found=1; next }
    found && /^\* /         { sub(/^\* /, ""); print; exit }
' "$RESULTS_FILE")

LATEST_PHD=$(printf '%s' "$LATEST_PHD" | sed -E '
    s/^([^,]+), \[(.*)\]\([^)]*\),[[:space:]]*\(\*\*[0-9]{4}\*\*\).*$/\1 — \2/
    s/^([^,]+), (.*),[[:space:]]*\(\*\*[0-9]{4}\*\*\).*$/\1 — \2/
')

echo "Most recent publication: $LATEST_PUB"
echo "Most recent PhD:         $LATEST_PHD"

if [[ -n "$LATEST_PUB" ]]; then
    sed -i "s|___XYZ_PUBLICATION___|$(prep "$LATEST_PUB")|" "$MAINPAGE"
else
    echo "WARNING: could not determine latest publication"
fi

if [[ -n "$LATEST_PHD" ]]; then
    sed -i "s|___XYZ_PHD___|$(prep "$LATEST_PHD")|" "$MAINPAGE"
else
    echo "WARNING: could not determine latest PhD"
fi
