# MaCh3

\htmlonly[block]
<div class="mach3-intro">

<div class="mach3-intro-text">
\endhtmlonly

## Introduction

The Markov Chain 3 flavour is a framework born in 2013 as a Bayesian MCMC fitter for [T2K](https://t2k-experiment.org/pl/) oscillation analysis.
It has now been used for multiple T2K Oscillation analyses both at the Near and Far detectors throughout the years and is also used by the DUNE and HK oscillation analysis groups as well as for joint fits between T2K and NOvA and T2K and SK's atmospheric data.

The framework has also evolved to allow non MCMC modules to interrogate the likelihoods implemented.

If something is unclear please contact us via
- [Slack](https://t2k-experiment.slack.com/archives/C06EM0C6D7W/p1705599931356889)
- [Discussions](https://github.com/mach3-software/MaCh3/discussions)
- [Indico](https://indico.global/category/1289/) If you need a password, please reach out to MaCh3-leadership for access.

\htmlonly
</div>

<div class="mach3-mcmc-demo">

<video controls autoplay loop muted playsinline>
    <source src="MCMC_Example.mp4" type="video/mp4">
</video>

<div class="mach3-mcmc-caption">
    Markov-Chain Monte Carlo sampling
</div>

</div>

</div>
\endhtmlonly


## Start Here

\htmlonly[block]

<div class="mach3-start-grid">

<a class="mach3-start-card" href="https://github.com/mach3-software/MaCh3Tutorial">
    <div class="mach3-start-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">🎓</span>
            Tutorial
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-start-card-text">
        Run a dummy experiment and quickly explore MaCh3's functionality.
    </div>
</a>

<a class="mach3-start-card" href="userguide.html">
    <div class="mach3-start-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">📖</span>
            User Guide
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-start-card-text">
        Learn the concepts behind Markov Chain Monte Carlo and Bayesian analysis.
    </div>
</a>

<a class="mach3-start-card" href="modules.html">
    <div class="mach3-start-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">🧩</span>
            API Reference
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-start-card-text">
        Understand the high-level structure of the MaCh3 framework.
    </div>
</a>

<a class="mach3-start-card" href="ReleaseNotes.html">
    <div class="mach3-start-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">📝</span>
            Release Notes
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-start-card-text">
        Find when features and changes were introduced.
    </div>
</a>

</div>

<div class="mach3-extra-grid">

<a class="mach3-extra-card" href="citelist.html">
    <div class="mach3-extra-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">📚</span>
            Bibliography
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-extra-card-text">
        Browse the scientific references cited throughout the MaCh3 documentation.
    </div>
</a>

<a class="mach3-extra-card" href="./pyMaCh3/mainpage.html">
    <div class="mach3-extra-card-title">
        <span class="mach3-card-heading">
            <span class="mach3-card-icon">🐍</span>
            Python Interface
        </span>
        <span class="mach3-card-arrow">→</span>
    </div>
    <div class="mach3-extra-card-text">
        Explore pyMaCh3 and learn how to interact with MaCh3 from Python.
    </div>
</a>

</div>

\endhtmlonly

## Recent results & publications

\htmlonly
<div class="mach3-recent-results">
    <table>
        <tr>
            <td class="mach3-recent-results-icon">
                <svg width="60" height="60" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
                    <path d="M14 2H6C4.9 2 4 2.9 4 4V20C4 21.1 4.9 22 6 22H18C19.1 22 20 21.1 20 20V8L14 2Z"
                          stroke="var(--mach3-blue)" stroke-width="2"
                          stroke-linecap="round" stroke-linejoin="round"/>
                    <path d="M14 2V8H20"
                          stroke="var(--mach3-blue)" stroke-width="2"
                          stroke-linecap="round" stroke-linejoin="round"/>
                    <path d="M16 13H8"
                          stroke="var(--mach3-blue)" stroke-width="2"
                          stroke-linecap="round" stroke-linejoin="round"/>
                    <path d="M16 17H8"
                          stroke="var(--mach3-blue)" stroke-width="2"
                          stroke-linecap="round" stroke-linejoin="round"/>
                    <path d="M10 9H8"
                          stroke="var(--mach3-blue)" stroke-width="2"
                          stroke-linecap="round" stroke-linejoin="round"/>
                </svg>
            </td>
            <td class="mach3-recent-results-text">
                <div>14+ papers using MaCh3</div>
                <div>34+ theses and dissertations</div>
                <div>2013 → present</div>
            </td>
            <td class="mach3-recent-results-link-icon">
                <svg width="70" height="30" viewBox="0 0 70 30" fill="none"
                    xmlns="http://www.w3.org/2000/svg">
                    <!-- Neutrino oscillation gradually becoming a direction -->
                    <path d="
                        M3 15
                        C7 15, 9 6, 13 6
                        C17 6, 19 24, 23 24
                        C27 24, 29 7, 33 7
                        C37 7, 39 23, 43 23
                        C47 23, 49 12, 53 12
                        C56 12, 58 15, 61 15
                    "
                    stroke="var(--mach3-blue)"
                    stroke-width="2"
                    stroke-linecap="round"
                    stroke-linejoin="round"
                    fill="none"/>
                    <!-- Arrowhead, continuing from the curve -->
                    <path d="M56 10L62 15L56 20"
                        stroke="var(--mach3-blue)"
                        stroke-width="2"
                        stroke-linecap="round"
                        stroke-linejoin="round"/>
                </svg>
            </td>
            <td class="mach3-recent-results-link">
                <a href="ResultsPublications.html">
                    See the full list of results and publications
                </a>
            </td>
        </tr>
    </table>
</div>
\endhtmlonly
