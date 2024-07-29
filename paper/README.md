# Intructions for the paper

NOTES:

- For the JOSS submission, the paper is written in markdown format and saved as paper.md. The bibliography is saved in paper.bib.

- For Arxiv submission, one can generate the paper in latex format using pandoc. The instructions are written in the Makefile and template of the paper given in latex.template file. Run `make` command in the terminal to generate the latex file and the corresponding pdf file. But the bibliography generated might not be compatible with the Arxiv submission. So, the bibliography should be manually adjusted and saved in the arxiv directory.

- Instruction for changing the bibliography in the latex file:

The following lines in the original latex file,

```
\section*{References}\label{references}
\addcontentsline{toc}{section}{References}

\phantomsection\label{refs}
\begin{CSLReferences}{1}{0}
\bibitem[\citeproctext]{ref-Abbott2021}
Abbott, R., T. D. Abbott, S. Abraham, F. Acernese, K. Ackley, A. Adams,
C. Adams, et al. 2021. {``Search for Lensing Signatures in the
Gravitational-Wave Observations from the First Half of LIGO--Virgo's
Third Observing Run.''} \emph{The Astrophysical Journal} 923 (1): 14.
\url{https://doi.org/10.3847/1538-4357/ac23db}.

\end{CSLReferences}
```

should be changed to,

```
\hypertarget{references}{%
\section*{References}\label{references}}
\addcontentsline{toc}{section}{References}

\hypertarget{refs}{}
\leavevmode\hypertarget{ref-Abbott2021}
Abbott, R., T. D. Abbott, S. Abraham, F. Acernese, K. Ackley, A. Adams,
C. Adams, et al. 2021. {``Search for Lensing Signatures in the
Gravitational-Wave Observations from the First Half of LIGO--Virgo's
Third Observing Run.''} \emph{The Astrophysical Journal} 923 (1): 14.
\url{https://doi.org/10.3847/1538-4357/ac23db}.
```