#
#
# This script is to build a qqplot first for the normal distribution, the idea come from UWA and the script ReliabilitySupportFms
#.which has problems with some data set
#
#

get_normal <-function(dataset){
        size_data <-length(dataset)
        na_vector <-is.na(dataset)
        if( sum(na_vector) > 0){
                has_na <-TRUE
                dataset <-dataset[!na_vector]
        }
        mean_data <-mean(dataset)





}



test_1 <-c(1.2,1.4,1.6,1.6,1.8, 1.8, 1.8, 2, 2.1)
