\ProvidesClass{tudothesis}[2015/08/31 v0.9 wrapper to scrbook for thesis at TUDo]
% Option for tu colors:
\newif\if@tucolor\@tucolorfalse
\DeclareOption{tucolor}{\@tucolortrue}

% All other Options go to KOMA-Script
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrbook}}
\PassOptionsToPackage{bibliography=totoc}{biblatex}
\ProcessOptions\relax

\LoadClass[
  headsepline,            % line  between header and text body
  bibliography=totoc,     % bibliography as unnumbered chapter in toc
  numbers=noenddot,       % no dot after figure/table number
  captions=tableheading,  % correct spacing for table headings
  titlepage=firstiscover, % symmetrical margins on titlepage
  headings=normal         % size of chapter headings slightly smaller
]{scrbook}


% Define our beloved tu colour:
\RequirePackage{xcolor}
\xdefinecolor{tugreen}{RGB}{132, 184, 25}
\colorlet{tulight}{tugreen!30!white}
\colorlet{tumed}{tugreen!90!black}
\colorlet{tudark}{tugreen!70!black}

\RequirePackage{caption}
\RequirePackage{subcaption}

% add some Color to structural elements if option tucolor was given
\RequirePackage{scrlayer-scrpage}
\pagestyle{scrheadings}

\KOMAoption{headsepline}{0.5pt}
\if@tucolor%
  \addtokomafont{title}{\color{tugreen}}          % title on titlepage
  \addtokomafont{chapter}{\color{tugreen}}     % headlines of the chapters, sections etc. 
  \setkomafont{headsepline}{\color{tugreen}}        % line between header and text body
  \DeclareCaptionFont{tugreen}{\color{tugreen}}
  \captionsetup{%
              labelfont={bf,tugreen},     % bold label in tu colour
              format=plain,               % no indentation under the label
              width=0.9\textwidth,        % slightly smaller width than text
              font=small,                 % font slightly smaller than text
            }
\else
  \captionsetup{%
              labelfont=bf,               % bold label
              format=plain,               % no indentation under the label
              width=0.9\textwidth,        % slightly smaller width than text
              font=small,                 % font slightly smaller than text
            }
\fi

% equation numbers of form chapter.equationnumber
\RequirePackage{amsmath}%
\renewcommand{\theequation}{\thechapter{}.\arabic{equation}}%
\numberwithin{equation}{chapter}%

% title page elements
\newcommand\thethesisclass{}
\newcommand\thetitle{}
\newcommand\thesubtitle{}
\newcommand\theauthor{}
\newcommand\thebirthplace{}
\newcommand\thedate{}
\newcommand\thechair{}
\newcommand\thedivision{}
\newcommand\theuniversity{}
\newcommand\thesisclass[1]{\renewcommand\thethesisclass{#1}}
\renewcommand\title[1]{\renewcommand\thetitle{#1}}
\renewcommand\subtitle[1]{\renewcommand\thesubtitle{#1}}
\renewcommand\author[1]{\renewcommand\theauthor{#1}}
\newcommand\birthplace[1]{\renewcommand\thebirthplace{geboren in #1}}
\renewcommand\date[1]{\renewcommand\thedate{#1}}
\newcommand\chair[1]{\renewcommand\thechair{#1}}
\newcommand\division[1]{\renewcommand\thedivision{#1}}
\newcommand\university[1]{\renewcommand\theuniversity{#1}}
\renewcommand{\subject}{Arbeit zur Erlangung des akademischen Grades \\ \thethesisclass}
\renewcommand{\publishers}{\thechair \\ \thedivision \\ \theuniversity}


\RequirePackage{eso-pic}
\RequirePackage{graphicx}

\newcommand{\TitlePageBackground}[1]{%
  \AddToShipoutPicture{%
    \includegraphics[width=\paperwidth,height=\paperheight,keepaspectratio]{#1}%
  }%
}

% custom title page
\renewcommand{\maketitle}{
  \begin{titlepage}
    \newgeometry{top=2.5cm, bottom=5cm, left=4.2cm, right=1.8cm}%
    \thispagestyle{empty}%
    \TitlePageBackground{Graphics/titlepage.pdf}%
    %\includegraphics[width=0.6\textwidth]{logos/tu-logo.pdf}%
    \begin{center}
      \vspace{2cm}%
      \rule{\textwidth}{1.2pt}\\
      \vspace{0.9cm}%
      \Huge\textsf{\textbf{\color{black}\thetitle}}\\
      \vspace{0.2cm}%
      \Huge{\textbf{--}}\\
      \vspace{0.2cm}%
      \Huge\textsf{\textbf{\color{black}\thesubtitle}}\\
      \vspace{0.5cm}%
      \rule{\textwidth}{1.2pt}\\
      \vspace{2.6cm}%
      \Large\theauthor\\
      \vspace{0.2cm}%
      \thedate\\
      \vfill%
    \end{center}
    \restoregeometry%
    \ClearShipoutPicture% remove background
  \end{titlepage}
}

% corrector page elements
\newcommand\thesubmissiondate{}
\newcommand\thefirstcorrector{}
\newcommand\thesecondcorrector{}
\newcommand\thethirdcorrector{}
\newcommand\firstcorrector[1]{\renewcommand\thefirstcorrector{#1}}
\newcommand\secondcorrector[1]{\renewcommand\thesecondcorrector{#1}}
\newcommand\thirdcorrector[1]{\renewcommand\thethirdcorrector{#1}}
\newcommand\submissiondate[1]{\renewcommand\thesubmissiondate{#1}}

% page for the correctors
\newcommand{\makecorrectorpage}{%
  \thispagestyle{empty}
  \vspace*{\fill}%
  \begin{tabular}{@{}l@{\hspace{2em}}l@{}}
      Erstgutachter:  & \thefirstcorrector\\
      Zweitgutachter: & \thesecondcorrector\\
      Drittgutachter: & \thethirdcorrector\\
      Abgabedatum:    & \thesubmissiondate%
  \end{tabular}
  \newpage
}

