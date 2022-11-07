# counts accessor for InteractionArray

    Code
      counts(iarr, TRUE)
    Output
      <11 x 11 x 10 x 2> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      
      ,,10,LEUK_HEK_PJA30_inter_30.hic
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
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      
      ,,3,LEUK_HEK_PJA30_inter_30.hic
               23760000 23765000 23770000 ... 23805000 23810000
      23700000        0        0        0   .        0        0
      23705000        0        0        0   .        0        0
      23710000        0        0        0   .        0        0
      23715000        0        0        0   .        0        0
      23720000        0        0        0   .        0        0
      23725000        0        0        0   .        0        0
      23730000        0        0        0   .        0        0
      23735000        0        0        0   .        0        0
      23740000        0        0        0   .        0        0
      23745000        0        0        0   .        0        0
      23750000        0        0        0   .        0        0

---

    Code
      counts(iarr[3:4, 1])
    Output
      <11 x 11 x 2 x 1> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
               23760000 23765000 23770000 ... 23805000 23810000
      23700000        0        0        0   .        0        0
      23705000        0        0        0   .        0        0
      23710000        0        0        0   .        0        0
      23715000        0        0        0   .        0        0
      23720000        0        0        0   .        0        0
      23725000        0        0        0   .        0        0
      23730000        0        0        0   .        0        0
      23735000        1        0        0   .        0        0
      23740000        0        0        0   .        1        1
      23745000        0        1        0   .        0        0
      23750000        1        0        0   .        0        0
      
      ...
      
      ,,2,LEUK_HEK_PJA27_inter_30.hic
                128645000 128650000 128655000 ... 128690000 128695000
      128140000         0         0         0   .         0         0
      128145000         0         0         0   .         0         0
      128150000         0         0         0   .         0         0
      128155000         0         0         0   .         0         0
      128160000         0         0         0   .         0         0
      128165000         0         0         0   .         0         0
      128170000         0         0         0   .         0         0
      128175000         0         0         0   .         0         0
      128180000         0         0         0   .         0         0
      128185000         0         0         0   .         0         0
      128190000         0         0         0   .         0         0

---

    Code
      counts(iarr[1:7, 1:2])
    Output
      <11 x 11 x 7 x 2> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      
      ,,7,LEUK_HEK_PJA30_inter_30.hic
                101510000 101515000 101520000 ... 101555000 101560000
      101055000         0         0         0   .         0         0
      101060000         0         0         0   .         0         0
      101065000         0         0         0   .         0         0
      101070000         0         0         0   .         0         0
      101075000         0         0         0   .         0         0
      101080000         0         0         0   .         0         0
      101085000         0         0         0   .         0         0
      101090000         0         0         0   .         0         0
      101095000         0         0         0   .         0         0
      101100000         0         0         0   .         0         0
      101105000         0         0         0   .         0         0

---

    Code
      counts(iarr[1:7, 1])
    Output
      <11 x 11 x 7 x 1> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      
      ,,7,LEUK_HEK_PJA27_inter_30.hic
                101510000 101515000 101520000 ... 101555000 101560000
      101055000         0         0         0   .         0         0
      101060000         0         0         0   .         0         0
      101065000         0         0         0   .         0         0
      101070000         0         0         0   .         0         0
      101075000         0         0         0   .         0         0
      101080000         0         0         0   .         0         0
      101085000         0         0         0   .         0         0
      101090000         0         0         0   .         0         0
      101095000         0         0         0   .         0         0
      101100000         0         0         0   .         0         0
      101105000         0         0         0   .         0         0

---

    Code
      counts(iarr[1, 1:2])
    Output
      <11 x 11 x 1 x 2> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      
      ,,1,LEUK_HEK_PJA30_inter_30.hic
               14740000 14745000 14750000 ... 14785000 14790000
      14435000        0        0        0   .        0        0
      14440000        0        0        0   .        0        0
      14445000        0        0        0   .        0        0
      14450000        0        0        0   .        0        0
      14455000        0        0        0   .        0        0
      14460000        0        0        0   .        0        1
      14465000        0        0        0   .        0        0
      14470000        0        0        0   .        0        0
      14475000        0        0        0   .        0        0
      14480000        0        0        0   .        0        0
      14485000        0        0        0   .        0        0

---

    Code
      counts(iarr[1, 1])
    Output
      <11 x 11 x 1 x 1> array of class DelayedArray and type "double":
      ,,1,LEUK_HEK_PJA27_inter_30.hic
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
      colnames(2): LEUK_HEK_PJA27_inter_30.hic LEUK_HEK_PJA30_inter_30.hic
      colData names(2): files fileNames
      type: GInteractions
      regions: 20

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
      colnames(2): LEUK_HEK_PJA27_inter_30.hic LEUK_HEK_PJA30_inter_30.hic
      colData names(2): files fileNames
      type: GInteractions
      regions: 20

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
      colnames(1): LEUK_HEK_PJA27_inter_30.hic
      colData names(2): files fileNames
      type: GInteractions
      regions: 20

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
      colnames(1): LEUK_HEK_PJA30_inter_30.hic
      colData names(2): files fileNames
      type: GInteractions
      regions: 20

