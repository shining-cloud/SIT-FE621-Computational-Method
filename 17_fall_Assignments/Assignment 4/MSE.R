D1 
data1


actual <- list(D1)
predicted <- list(data1)

for(i in length(actual)){
    cat("MSE ")
    cat(mean((actual[[i]] - predicted[[i]])^2))
}
    mean((D1-data1)^2)
