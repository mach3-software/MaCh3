#dnf install -y texlive-scheme-basic texlive-latex
#dnf install -y mathjax
#dnf install -y perl
#dnf install -y doxygen
git clone https://github.com/mach3-software/MaCh3.wiki.git Doc/User-Guide
BASE_URL="https://raw.githubusercontent.com/KSkwarczynski/Python/main/MaCh3/animation/"
curl -L "$BASE_URL/MCMC_Example.mp4" -o Doc/MaCh3Web/MCMC_Example.mp4
curl -L "$BASE_URL/posterior_predictive_Example.mp4" -o Doc/MaCh3Web/posterior_predictive_Example.mp4
curl -L "$BASE_URL/umbrella_sampling_delta_cp.mp4" -o Doc/MaCh3Web/umbrella_sampling_delta_cp.mp4
