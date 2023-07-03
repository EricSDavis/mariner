# Interpolation works with single row/column

    Code
      counts(regularize(ija, ndim = c(3, 3), scale = FALSE))
    Message <simpleMessage>
      / Reading and realizing block 1/3 of file 1/1 ... 
      OK
      \ Processing it ... 
      OK
      / Reading and realizing block 2/3 of file 1/1 ... 
      OK
      \ Processing it ... 
      OK
      / Reading and realizing block 3/3 of file 1/1 ... 
      OK
      \ Processing it ... 
      OK
    Output
      <3 x 3 x 3 x 1> DelayedArray object of type "double":
      ,,1,FS
           [,1] [,2] [,3]
      [1,]   53   15    5
      [2,]   53   15    5
      [3,]   53   15    5
      
      ,,2,FS
           [,1] [,2] [,3]
      [1,]   53   53   53
      [2,]   15   15   15
      [3,]    5    5    5
      
      ,,3,FS
           [,1] [,2] [,3]
      [1,]   53   53   53
      [2,]   53   53   53
      [3,]   53   53   53
      

---

    Code
      counts(regularize(ija, ndim = c(1, 3), scale = FALSE, verbose = FALSE))
    Output
      <1 x 3 x 3 x 1> DelayedArray object of type "double":
      ,,1,FS
           [,1] [,2] [,3]
      [1,]   53   15    5
      
      ,,2,FS
           [,1] [,2] [,3]
      [1,]   53   53   53
      
      ,,3,FS
           [,1] [,2] [,3]
      [1,]   53   53   53
      

---

    Code
      counts(regularize(ija, ndim = c(3, 1), scale = FALSE, verbose = FALSE))
    Output
      <3 x 1 x 3 x 1> DelayedArray object of type "double":
      ,,1,FS
           [,1]
      [1,]   53
      [2,]   53
      [3,]   53
      
      ,,2,FS
           [,1]
      [1,]   53
      [2,]   15
      [3,]    5
      
      ,,3,FS
           [,1]
      [1,]   53
      [2,]   53
      [3,]   53
      

