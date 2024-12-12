%% myhomework.cls
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{math_hw}[2024/02/04 Custom Homework Class]

\LoadClass[12pt, letterpaper]{article}

\RequirePackage[T1]{fontenc}
\RequirePackage{tgpagella}
\RequirePackage{mathtools}
\RequirePackage{geometry}
\RequirePackage{microtype}
\RequirePackage{titlesec}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{setspace}
\RequirePackage{hyperref}
\RequirePackage{fancyhdr}
\RequirePackage{graphicx}
\RequirePackage{float}
\RequirePackage{xcolor}
\definecolor{LightGray}{gray}{0.95}
\RequirePackage{apptools}
\RequirePackage{minted}
\RequirePackage{authblk}
\RequirePackage{units}
% \RequirePackage[inkscapeformat=png]{svg}
\RequirePackage[labelfont=bf]{caption}
\RequirePackage[noabbrev,capitalize]{cleveref}
\RequirePackage{amsthm}
\RequirePackage{bm}
\RequirePackage{subcaption}

\renewcommand\qedsymbol{\blacksquare}

\geometry{margin=1in}

\titleformat{\section}{\normalfont\Large\itshape}{\textbf{Part~\thesection: }}{1em}{}
\titleformat{\subsection}{\normalfont\large}{\textbf{\thesubsection: }}{1em}{}
\AtAppendix{\titleformat{\section}{\normalfont\Large\itshape}{\textbf{Appendix~\thesection: }}{1em}{}}

\hypersetup{
    colorlinks=true,
    linkcolor=black,
    filecolor=magenta,
    urlcolor=blue,
    citecolor=black
}

\renewcommand{\topfraction}{0.85}
\renewcommand{\textfraction}{0.1}
\renewcommand{\floatpagefraction}{0.7}
\floatplacement{figure}{H}

\renewcommand{\headrulewidth}{0pt}
\pagestyle{fancy}
\fancyhf{}
\setlength{\headheight}{20pt}

\onehalfspacing

\renewcommand{\maketitle}{
  \begin{titlepage}
    \begin{center}
      {\LARGE{\@title}}\par
      \vspace{2em}

      \vspace{8em}
      \lineskip .5em
      {\large
        \begin{tabular}[t]{c}
          \@author
        \end{tabular}\par
      }
        \vfill
        {\large\textbf{\@date}}\par
    \end{center}
  \end{titlepage}
}

\endinput
