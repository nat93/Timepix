TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += "/home/andrii/root/include"
#INCLUDEPATH += "/home/anatochi/root-6.10.02/include"

LIBS += -L/home/andrii/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread
#LIBS += -L/home/anatochi/root-6.10.02/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread



HEADERS += \
    .gitignore

DISTFILES += \
    run.bash \
    autorun_1.bash \
    autorun_2.bash \
    README.md

SOURCES += \
    analysis_trk.cc \
    trackreco_trk.cc \
    analysis_common.cc \
    ascii2root_common.cc \
    convert_common.cc \
    csv2root_common.cc \
    ch_eff/sps_md_analysis_2.C \
    ch_eff/sps_md_analysis_1.C \
    ch_eff/sps_md_analysis_3.C \
    ch_eff/sps_md_analysis_4.C \
    ch_eff/sps_md_analysis_1.C \
    ch_eff/sps_md_analysis_1_double_channeling.C \
    ch_eff/sps_md_analysis_all_double_channeling.C \
    histos_may.C \
    histos_september.C \
    ch_eff/test_05_09_2018.C \
    ch_eff/sps_md_analysis_17_09_2018.C \
    ch_eff/sps_md_analysis_17_09_2018_rp0_issue.C \
    ch_eff/sps_md_analysis_17_09_2018_double_channeling.C \
    plot_tree.C \
    ch_eff/sps_md_analysis_07_11_2018_double_channeling.C \
    ch_eff/sps_md_analysis_17_10_2017_thesis.cc \
    ch_eff/sps_md_analysis_17_09_2018_thesis.cc \
    ch_eff/sps_md_analysis_07_11_2018_thesis.cc \
    ch_eff/ion_analysis.cc \
    ch_eff/plot_taratin_extraction.C \
    ch_eff/plot_ion_cluster.C



