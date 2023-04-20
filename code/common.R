#
#
# Commen for the MTBF 
#
#
library(comparator)

#
#
# distance is the calculation for token text list input. The distance is Monge Elkan for the moment
#
# input:  a vector of list where each list is a set of token 
#
# ouptut: squre matrix with the Monge Elkan distance 
#
distance <-function( vector, distance= NULL){
  assert::assert(class(vector) == "list")
  if( is.null(distance)){
    distance <- "Monge-Elkan"
  }
  # We build the matrix
  total_y <-length(vector)
  matrix_distance <-matrix( ncol=total_y)
  for(index_y in 1:total_y){
    x_calc <-numeric()
    for(index_x in 1:total_y){
      x_calc[index_x] <-MongeElkan(agg_function= hmean)(list(vector[[index_y]]), list(vector[[index_x]]))
    }
    matrix_distance <-rbind(matrix_distance, x_calc)
  }
  return(matrix_distance)
}
  