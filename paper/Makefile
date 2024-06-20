default: pdf

MDFILE=paper.md
BIBFILE=paper.bib
LATEXFILE=paper.tex

# some variables for the paper
REPOSITORY="https://github.com/hemantaph/ler"
REVIEW=""
ARCHIVEDOI="https://zenodo.org/badge/latestdoi/626733473"
FORMATTEDDOI="10.xxxxxx/draft"
YEAR="2024"
TITLE="$ler$ : LVK (LIGO-Virgo-KAGRA collaboration) Event (compact-binary mergers) Rate calculator and simulator"
AUTHOR="Phurailatpam Hemantakumar, Anupreeta More, Harsh Narola, Ng Chung Yin (Leo), Justin Janquart, Chris Van Den Broeck, Otto Akseli HANNUKSELA, Neha Singh, and David Keitel"
ISSUE=""
VOLUME=""
PAGE=""
JOURNALNAME="Journal of Open Source Software"
PUBLISHED="19 June 2024"
SUBMITTED="19 June 2024"
PAPERURL=""
EDITOR=""
EDITORURL=""
REVIEWER1=""
REVIEWER2=""

# get JOSS template and logo
latex.template logo.png:
	curl "https://raw.githubusercontent.com/openjournals/whedon/master/resources/joss/$@" > $@

# PREPRINT=true
ifdef PREPRINT
  LOGO=logo-empty.png
  SUBMITBIT="Submitted to The "
else
  LOGO=logo.png
  SUBMITBIT=""
endif


# use pandoc with native citation processing to create LaTeX source
latex: latex.template $(LOGO) $(MDFILE) $(BIBFILE)
	pandoc --citeproc --bibliography=$(BIBFILE) -V repository=$(REPOSITORY) -V journal_name=$(SUBMITBIT)$(JOURNALNAME) -V paper_url=$(PAPERURL) -V archive_doi=$(ARCHIVEDOI) -V review_issue_url=$(REVIEW) -V formatted_doi=$(FORMATTEDDOI) -V logo_path=$(LOGO) -V year=$(YEAR) -V paper_title=$(TITLE) -V citation_author=$(AUTHOR) -V issue=$(ISSUE) -V volume=$(VOLUME) -V page=$(PAGE) -V submitted=$(SUBMITTED) -V published=$(PUBLISHED) -V editor_name=$(EDITOR) -V editor_url=$(EDITORURL) -V reviewers=$(REVIEWER1) -V reviewers=$(REVIEWER2) -V graphics=true -V papersize=a4paper -V geometry:margin=1in -t latex -s $(MDFILE) -o $(LATEXFILE) --template latex.template



pdf: latex
	pdflatex paper.tex
	pdflatex paper.tex

clean:
	rm paper.bcf paper.aux paper.log paper.out paper.run.xml

clean_downloads:
	rm latex.template logo.png