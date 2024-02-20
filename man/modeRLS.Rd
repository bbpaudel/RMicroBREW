\name{modeRLS}
\alias{modeRLS}
\title{Quantify modes of different budding}
\usage{
modeRLS(data, mode_Frames, plot_it, save_file)
}
\description{
This extracts data required for RLS of yeast cells, and performs statistical analysis to identify distinct modes of budding. The user can define how many modes they expect to see in the data.
}
\examples{
mode_output = modeRLS(data = output[[2]], mode_Frames = 2, plot_it = "yes", save_file = "file_name")
}
\arguments{
\item{data}{output[[2]] table from MicroBrew_count function, which contains information as the peak_location, which shape it came from, mother cell entry and death if identified.}
\item{mode_Frames}{number of modes expected in the data}
\item{plot_it}{either "yes" or "no" from the analysis, this will plot a boxplot with the peak_height in y-axis and classification or modes in the x-axis}
\item{save_file}{whether to save the file or the results of the analysis, and what name to use}
}
\value{
The function uses Mclust, R package to fit mixed gaussian distributions with the peak height as variables. It outputs the proportions of data points in distinct modes, and also saves the classification variable back to the original data to be used for other detailed analysis.
}
