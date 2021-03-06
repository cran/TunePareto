\name{NEWS}
\title{NEWS file for the \pkg{TunePareto} package}

\section{Changes in version 2.5}{

\subsection{Modified functions export}{

\itemize{

\item{Fix in exporting predict.TuneParetoModel.}

}
}
}

\section{Changes in version 2.4}{

\subsection{Modified functions}{

\itemize{

\item{Fix in predict.TuneParetoModel preventing returning empty prediction vectors.}

}
}
}

\section{Changes in version 2.3}{

\subsection{New functions}{

\itemize{

\item{New objective functions: 
\code{reclassAccuracy},
\code{reclassRecall},
\code{reclassPrecision},
\code{reclassFallout},
\code{reclassFalsePositive},
\code{reclassMiss},
\code{reclassFalseNegative},
\code{reclassTruePositive},
\code{reclassTrueNegative},
\code{reclassPPV},
\code{reclassNPV},
\code{cvAccuracy},
\code{cvRecall},
\code{cvPrecision},
\code{cvFallout},
\code{cvFalsePositive},
\code{cvMiss},
\code{cvFalseNegative},
\code{cvTruePositive},
\code{cvTrueNegative},
\code{cvPPV},
\code{cvNPV}}
}

}
\subsection{Modified functions}{

\itemize{

\item{All objective functions now have a parameter \code{saveModel} that stores the classification model in the result object.}

}
}
}
