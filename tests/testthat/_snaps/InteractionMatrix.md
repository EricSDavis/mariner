# show method for InteractionArray

    Code
      show(imat)
    Output
      class: InteractionMatrix 
      dim: count matrix with 3 interactions and 2 file(s)
      metadata(3): binSize norm matrix
      assays(1): counts
      rownames: NULL
      rowData names(0):
      colnames(2): FS WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 4

---

    Code
      show(imat[1:2, ])
    Output
      class: InteractionMatrix 
      dim: count matrix with 2 interactions and 2 file(s)
      metadata(3): binSize norm matrix
      assays(1): counts
      rownames: NULL
      rowData names(0):
      colnames(2): FS WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 4

---

    Code
      show(imat[1:2, 1])
    Output
      class: InteractionMatrix 
      dim: count matrix with 2 interactions and 1 file(s)
      metadata(3): binSize norm matrix
      assays(1): counts
      rownames: NULL
      rowData names(0):
      colnames(1): FS
      colData names(2): files fileNames
      type: GInteractions
      regions: 4

---

    Code
      show(imat[, 2])
    Output
      class: InteractionMatrix 
      dim: count matrix with 3 interactions and 1 file(s)
      metadata(3): binSize norm matrix
      assays(1): counts
      rownames: NULL
      rowData names(0):
      colnames(1): WT
      colData names(2): files fileNames
      type: GInteractions
      regions: 4

