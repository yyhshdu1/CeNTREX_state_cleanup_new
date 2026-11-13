# LaTeX environment for editing the state-cleanup report (docs/state_cleanup/main.tex) in GitHub Codespaces
# or a local dev container. Only the TeX packages the report needs, not the full TeX Live (~5 GB).
FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        latexmk \
        lmodern \
        texlive-latex-base \
        texlive-latex-recommended \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-science \
        chktex \
        poppler-utils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
