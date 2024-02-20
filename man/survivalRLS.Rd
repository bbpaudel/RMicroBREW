\name{survivalRLS}
\alias{survivalRLS}
\title{Get data for survival analysis}
\usage{
survivalRLS(data_rls, numFrames, group_name, plot_it, save_file, color1)
}
\description{
This extracts data for Kaplan-Meier type plots. This function takes the output[[2]] of the MicroBrew_count function and estimates survival probability and returns a table with total budding events, and survival status.
}
\examples{
surv_data = survivalRLS(data_rls = data, numFrames = 526, group_name = "group_name", plot_it = "yes", save_file = "file_name", color1 = 'red')
}
\arguments{
\item{data_rls}{output[[2]] table from MicroBrew_count function, which contains information as the peak_location, which shape it came from, mother cell entry and death if identified.}
\item{numFrames}{number of frames to process, this number is critical to decide the survival or death status, if no events are identified before the numFrames, then it is assigned no event.}
\item{group_name}{name of the experimental group or the name of the images processed.}
\item{plot_it}{either "yes" or "no" to get the survival plot from the analysis.}
\item{save_file}{whether to save the results, provide the name for the file, "name_file"}
\item{color1}{what color to be used for survival plot}
}
\value{
The function outputs data extracted from MicroBrew_count that are used for survival plots. The output table contains information about the trap_id, cell entry, budding events from shape1, shape2, shape3, total budding events, death frame, live/dead status (0 is no event or live, 1 is event or censored or death), and group_id. The function also saves RLS_count data into the folder called MicroBREW_results.
}
