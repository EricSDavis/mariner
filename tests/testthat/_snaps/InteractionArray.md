# counts accessor for InteractionArray

    Code
      counts(iarr, TRUE)
    Output
      <11 x 11 x 10 x 2> array of class DelayedArray and type "double":
      ,,1,FS
               14740000 14745000 14750000 ... 14785000 14790000
      14435000        0        0        0   .        0        0
      14440000        0        0        0   .        0        0
      14445000        0        0        0   .        0        0
      14450000        0        0        0   .        0        0
      14455000        0        0        0   .        0        0
      14460000        0        1        0   .        0        0
      14465000        0        0        0   .        0        0
      14470000        0        0        0   .        0        0
      14475000        0        0        0   .        0        0
      14480000        0        0        0   .        0        0
      14485000        0        0        0   .        0        0
      
      ...
      
      ,,10,WT
                104045000 104050000 104055000 ... 104090000 104095000
      103335000         0         0         0   .         0         0
      103340000         0         0         0   .         0         0
      103345000         0         0         0   .         0         0
      103350000         0         0         0   .         0         0
      103355000         0         0         0   .         0         0
      103360000         0         0         0   .         0         0
      103365000         0         0         0   .         0         0
      103370000         0         0         0   .         0         0
      103375000         0         0         0   .         0         0
      103380000         0         0         0   .         0         0
      103385000         0         0         0   .         0         0

---

    Code
      counts(iarr[1:3, 1:2])
    Output
      <11 x 11 x 3 x 2> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      
      ...
      
      ,,3,WT
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

---

    Code
      counts(iarr[3:4, 1])
    Output
      <11 x 11 x 2 x 1> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     1     0   .     0     0
      [11,]     1     0     0   .     0     0
      
      ,,2,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

---

    Code
      counts(iarr[1:7, 1:2])
    Output
      <11 x 11 x 7 x 2> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      
      ...
      
      ,,7,WT
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

---

    Code
      counts(iarr[1:7, 1])
    Output
      <11 x 11 x 7 x 1> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      
      ...
      
      ,,7,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

---

    Code
      counts(iarr[1, 1:2])
    Output
      <11 x 11 x 1 x 2> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      
      ,,1,WT
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

---

    Code
      counts(iarr[1, 1])
    Output
      <11 x 11 x 1 x 1> array of class DelayedArray and type "double":
      ,,1,FS
             [,1]  [,2]  [,3] ... [,10] [,11]
       [1,]     0     0     0   .     0     0
       [2,]     0     0     0   .     0     0
        ...     .     .     .   .     .     .
      [10,]     0     0     0   .     0     0
      [11,]     0     0     0   .     0     0
      

# show method for InteractionArray

    Code
      show(iarr)
    Output
      class: InteractionArray 
      dim: 10 interaction(s), 2 file(s), 11x11 count matrix(es)
      metadata(3): binSize norm matrix
      assays(3): counts rownames colnames
      rownames: NULL
      rowData names(0):
      colnames(2): FS WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 20390

---

    Code
      show(iarr[1:3, ])
    Output
      class: InteractionArray 
      dim: 3 interaction(s), 2 file(s), 11x11 count matrix(es)
      metadata(3): binSize norm matrix
      assays(3): counts rownames colnames
      rownames: NULL
      rowData names(0):
      colnames(2): FS WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 20390

---

    Code
      show(iarr[1:3, 1])
    Output
      class: InteractionArray 
      dim: 3 interaction(s), 1 file(s), 11x11 count matrix(es)
      metadata(3): binSize norm matrix
      assays(3): counts rownames colnames
      rownames: NULL
      rowData names(0):
      colnames(1): FS
      colData names(2): files fileNames
      type: GInteractions
      regions: 20390

---

    Code
      show(iarr[, 2])
    Output
      class: InteractionArray 
      dim: 10 interaction(s), 1 file(s), 11x11 count matrix(es)
      metadata(3): binSize norm matrix
      assays(3): counts rownames colnames
      rownames: NULL
      rowData names(0):
      colnames(1): WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 20390

