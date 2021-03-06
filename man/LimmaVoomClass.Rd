% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voom.R
\name{LimmaVoomClass}
\alias{LimmaVoomClass}
\title{R6 class for limma voom analysis}
\value{

}
\description{
R6 class for limma voom analysis

R6 class for limma voom analysis
}
\examples{
result <- LimmaVoomClass$new(dge, groupcol, outdir)
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{LimmaVoomClass$new()}}
\item \href{#method-getResults}{\code{LimmaVoomClass$getResults()}}
\item \href{#method-clone}{\code{LimmaVoomClass$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LimmaVoomClass$new(dge, samples, outdir)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{outdir}}{directory to write DGE results to}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getResults"></a>}}
\if{latex}{\out{\hypertarget{method-getResults}{}}}
\subsection{Method \code{getResults()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LimmaVoomClass$getResults(p = NULL, padj = 0.05, logfc = 1.5)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{LimmaVoomClass$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
