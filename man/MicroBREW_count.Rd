\name{MicroBREW_count}
\alias{Function to extract information from MicroBREW analysis}
\title{provide summary of the function}
\usage{
MicroBrew_count(filename, save_file)
}
\description{
This function gets the count and extracts relevant information from output of accompanying MATLAB function processing the images acquired with Microfluidic chips. It takes in csv files as an input.
}
\examples{
results = MicroBREW_count("filename.csv", save_file = "result_name")
}
\value{
The function results to three table output. If output is stored as results, results[[1]] will contain information containing the following information: trap_id, cell_entry (to denote the frame mother cell enters the chip), number of events from shape1, shape2, and shape3 recorded separately.

The output, results[[2]] will contain information with identified peaks and their features, such as peak height (in pixels), the frames the peak was identified (peak_location), the frames the peak began (peak_begin), the frames the peak end (peak_end), which shape the peaks are associated with (which_shape), the frame of mother cell entry (mother_cell_entry), the trap_id, and the frame the mother cell is identified to have died if it died. If the mother cell does not die until the end of the experiment, the max number of frame will be recorded.

The output, results[[3]] will contain all the columns that the results from MATLAB contained along with new trap id, smoothened center intensity values, and smootheed values for three different shapes. It also contains information on the lag difference of the pixels between the frame with its earlier frame, which is used to detect when the mother cell enters the trap, and when it dies. The sum_value denotes the sum of the center intensity and the lag difference, used to identify the location of mother cell death. Column value_pass denotes if the algorithm identifies the mother cell entry or not.
}
