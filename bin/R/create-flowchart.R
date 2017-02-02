# library(DiagrammeR)
# DiagrammeR("
#                 graph TB
#            A(GC/ GC-MS)-->B(Peak Detection)
#            B-->C((GCalingR))
#            C-->D(1. Linear Alignment)
#            C-->E(2. Piecewise Alignment)
#            C-->F(3. Grouping)
#            D-->G(Alignment)
#            E-->G
#            F-->G
#            G-->H(Analysis)
#            style C fill:#1B9E77;
#            style D fill:#D95F02;
#            style E fill:#D95F02;
#            style F fill:#D95F02;
#            style G fill:#7570B3;
#            ",
#            height = 900, width = 600)

library(DiagrammeR)
a <- DiagrammeR(
    "graph TB
    A(GC/ GC-MS)-->B(Peak Detection)
    B-->C((GCalingR))
    C-->D(check_input)
    D-->E
    C-->E(align_chromatograms)
    E-->F(aligned data)
    F-->G(gc_heatmap)
    F-->H(plot)
    F-->I(print)
    style C fill:#1B9E77;
    style D fill:#D95F02;
    style E fill:#D95F02;
style G fill:#D95F02;
style H fill:#D95F02;
style I fill:#D95F02
    "
)
a
b <- DiagrammeR("graph LR
                A(aligned data)-->B(norm_peaks)
                style B fill:#D95F02")
c <- DiagrammeR("graph BT
                A(norm_peaks)-->B(vegan)
                B-->C(metaMDS)
                B-->D(adonis)
                B-->E(betadisper)
                style B fill:#1B9E77;
    style C fill:#D95F02;
                style D fill:#D95F02;
style E fill:#D95F02
")
