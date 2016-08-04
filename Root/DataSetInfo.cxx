/* Copyright 2016 Oleksandr Viazlo */

#include <MyAnalysis/DataSetInfo.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(DataSetInfo)

DataSetInfo::DataSetInfo(){
   
/// nEvents
m_nEventsMap[410016] = 997600;  /// PowhegPythiaEvtGen_P2012_Wt_dilepton_antitop
m_nEventsMap[410015] = 997200;	/// PowhegPythiaEvtGen_P2012_Wt_dilepton_top
m_nEventsMap[410014] = 4985600;	/// PowhegPythiaEvtGen_P2012_Wt_inclusive_antitop
m_nEventsMap[410013] = 4985800;	/// PowhegPythiaEvtGen_P2012_Wt_inclusive_top
m_nEventsMap[410012] = 4989800;	/// PowhegPythiaEvtGen_P2012_singletop_tchan_lept_antitop
m_nEventsMap[410011] = 4988200;	/// PowhegPythiaEvtGen_P2012_singletop_tchan_lept_top
m_nEventsMap[410000] = 19972600;/// PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
m_nEventsMap[361108] = 4541600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau
m_nEventsMap[361107] = 19964000;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu
m_nEventsMap[361106] = 19920600;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee
m_nEventsMap[361104] = 19970400;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusmunu
m_nEventsMap[361103] = 19922600;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu
m_nEventsMap[361101] = 19997800;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu
m_nEventsMap[361100] = 29898800;/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu
m_nEventsMap[361096] = 487000;	/// Sherpa_CT10_ZqqZll_SHv21_improved
m_nEventsMap[361095] = 973000;	/// Sherpa_CT10_WqqZvv_SHv21_improved
m_nEventsMap[361094] = 490000;	/// Sherpa_CT10_WqqZll_SHv21_improved
m_nEventsMap[361093] = 1952000;	/// Sherpa_CT10_WlvZqq_SHv21_improved
m_nEventsMap[361092] = 1983000;	/// Sherpa_CT10_WpqqWmlv_SHv21_improved
m_nEventsMap[361091] = 1977000;	/// Sherpa_CT10_WplvWmqq_SHv21_improved

/// FIXME in AMI nEvents equal to 366025 while in support note - 2000000
/// let's use value from support note for testing purpise
// m_nEventsMap[361088] = 366025;	/// Sherpa_CT10_lvvv
m_nEventsMap[361088] = 2000000;	/// Sherpa_CT10_lvvv

m_nEventsMap[361086] = 99400;	/// Sherpa_CT10_ZqqZll
m_nEventsMap[361084] = 249400;	/// Sherpa_CT10_WqqZll
m_nEventsMap[361068] = 5929600;	/// Sherpa_CT10_llvv
m_nEventsMap[361067] = 1197000;	/// Sherpa_CT10_lllvOFPlus
m_nEventsMap[361066] = 598600;	/// Sherpa_CT10_lllvSFPlus
m_nEventsMap[361065] = 898000;	/// Sherpa_CT10_lllvOFMinus
m_nEventsMap[361064] = 448800;	/// Sherpa_CT10_lllvSFMinus
m_nEventsMap[361063] = 990000;	/// Sherpa_CT10_llll
m_nEventsMap[301534] = 878400;	/// Pythia8EvtGen_A14NNPDF23LO_Wprime_munu_Flat_13TeV
m_nEventsMap[301533] = 997800;	/// Pythia8EvtGen_A14NNPDF23LO_Wprime_enu_Flat_13TeV
m_nEventsMap[301138] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_5000M
m_nEventsMap[301137] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_4500M5000
m_nEventsMap[301136] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_4000M4500
m_nEventsMap[301135] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_3500M4000
m_nEventsMap[301134] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_3000M3500
m_nEventsMap[301133] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2750M3000
m_nEventsMap[301132] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2500M2750
m_nEventsMap[301131] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2250M2500
m_nEventsMap[301130] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2000M2250
m_nEventsMap[301129] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1750M2000
m_nEventsMap[301128] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1500M1750
m_nEventsMap[301127] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1250M1500
m_nEventsMap[301126] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1000M1250
m_nEventsMap[301125] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_800M1000
m_nEventsMap[301124] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_600M800
m_nEventsMap[301123] = 99800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_400M600
m_nEventsMap[301122] = 149600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_250M400
m_nEventsMap[301121] = 249600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_180M250
m_nEventsMap[301120] = 499800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_120M180
m_nEventsMap[301118] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_5000M
m_nEventsMap[301117] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_4500M5000
m_nEventsMap[301116] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_4000M4500
m_nEventsMap[301115] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_3500M4000
m_nEventsMap[301114] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_3000M3500
m_nEventsMap[301113] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2750M3000
m_nEventsMap[301112] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2500M2750
m_nEventsMap[301111] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2250M2500
m_nEventsMap[301110] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2000M2250
m_nEventsMap[301109] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1750M2000
m_nEventsMap[301108] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1500M1750
m_nEventsMap[301107] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1250M1500
m_nEventsMap[301106] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1000M1250
m_nEventsMap[301105] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_800M1000
m_nEventsMap[301104] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_600M800
m_nEventsMap[301103] = 100000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_400M600
m_nEventsMap[301102] = 150000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_250M400
m_nEventsMap[301101] = 249400;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_180M250
m_nEventsMap[301100] = 499400;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_120M180
m_nEventsMap[301098] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_5000M
m_nEventsMap[301097] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_4500M5000
m_nEventsMap[301096] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_4000M4500
m_nEventsMap[301095] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_3500M4000
m_nEventsMap[301094] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_3000M3500
m_nEventsMap[301093] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2750M3000
m_nEventsMap[301092] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2500M2750
m_nEventsMap[301091] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2250M2500
m_nEventsMap[301090] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2000M2250
m_nEventsMap[301089] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1750M2000
m_nEventsMap[301088] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1500M1750
m_nEventsMap[301087] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1250M1500
m_nEventsMap[301086] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1000M1250
m_nEventsMap[301085] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_800M1000
m_nEventsMap[301084] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_600M800
m_nEventsMap[301083] = 100000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_400M600
m_nEventsMap[301082] = 149400;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_250M400
m_nEventsMap[301081] = 249800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_180M250
m_nEventsMap[301080] = 499600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_120M180
m_nEventsMap[301078] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_5000M
m_nEventsMap[301077] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_4500M5000
m_nEventsMap[301076] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_4000M4500
m_nEventsMap[301075] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_3500M4000
m_nEventsMap[301074] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_3000M3500
m_nEventsMap[301073] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2750M3000
m_nEventsMap[301072] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2500M2750
m_nEventsMap[301071] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2250M2500
m_nEventsMap[301070] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2000M2250
m_nEventsMap[301069] = 40000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1750M2000
m_nEventsMap[301068] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1500M1750
m_nEventsMap[301067] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1250M1500
m_nEventsMap[301066] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1000M1250
m_nEventsMap[301065] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_800M1000
m_nEventsMap[301064] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_600M800
m_nEventsMap[301063] = 99400;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_400M600
m_nEventsMap[301062] = 139600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_250M400
m_nEventsMap[301061] = 249600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_180M250
m_nEventsMap[301060] = 498800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_120M180
m_nEventsMap[301038] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_5000M
m_nEventsMap[301037] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_4500M5000
m_nEventsMap[301036] = 49000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_4000M4500
m_nEventsMap[301035] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_3500M4000
m_nEventsMap[301034] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_3000M3500
m_nEventsMap[301033] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2750M3000
m_nEventsMap[301032] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2500M2750
m_nEventsMap[301031] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2250M2500
m_nEventsMap[301030] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2000M2250
m_nEventsMap[301029] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1750M2000
m_nEventsMap[301028] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1500M1750
m_nEventsMap[301027] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1250M1500
m_nEventsMap[301026] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1000M1250
m_nEventsMap[301025] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_800M1000
m_nEventsMap[301024] = 11000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_600M800
m_nEventsMap[301023] = 99000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_400M600
m_nEventsMap[301022] = 149200;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_250M400
m_nEventsMap[301021] = 249600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_180M250
m_nEventsMap[301020] = 498600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_120M180
m_nEventsMap[301018] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_5000M
m_nEventsMap[301017] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_4500M5000
m_nEventsMap[301016] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_4000M4500
m_nEventsMap[301015] = 49600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_3500M4000
m_nEventsMap[301014] = 48600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_3000M3500
m_nEventsMap[301013] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2750M3000
m_nEventsMap[301012] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2500M2750
m_nEventsMap[301011] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2250M2500
m_nEventsMap[301010] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2000M2250
m_nEventsMap[301009] = 3600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1750M2000
m_nEventsMap[301008] = 49800;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1500M1750
m_nEventsMap[301007] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1250M1500
m_nEventsMap[301006] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1000M1250
m_nEventsMap[301005] = 50000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_800M1000
m_nEventsMap[301004] = 145000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_600M800
m_nEventsMap[301003] = 100000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_400M600
m_nEventsMap[301002] = 149400;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_250M400
m_nEventsMap[301001] = 249600;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_180M250
m_nEventsMap[301000] = 498000;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_120M180

/// xSec, [nb]
m_crossSectionMap[410016] = .0035814;	/// PowhegPythiaEvtGen_P2012_Wt_dilepton_antitop
m_crossSectionMap[410015] = .0035835;	/// PowhegPythiaEvtGen_P2012_Wt_dilepton_top
m_crossSectionMap[410014] = .033989;	/// PowhegPythiaEvtGen_P2012_Wt_inclusive_antitop
m_crossSectionMap[410013] = .034009;	/// PowhegPythiaEvtGen_P2012_Wt_inclusive_top
m_crossSectionMap[410012] = .025778;	/// PowhegPythiaEvtGen_P2012_singletop_tchan_lept_antitop
m_crossSectionMap[410011] = .043739;	/// PowhegPythiaEvtGen_P2012_singletop_tchan_lept_top
m_crossSectionMap[410000] = .69611;	/// PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_nonallhad
m_crossSectionMap[361108] = 1.9012;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau
m_crossSectionMap[361107] = 1.9012;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu
m_crossSectionMap[361106] = 1.9012;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee
m_crossSectionMap[361104] = 8.2826;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusmunu
m_crossSectionMap[361103] = 8.2826;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu
m_crossSectionMap[361101] = 11.306;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu
m_crossSectionMap[361100] = 11.306;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu
m_crossSectionMap[361096] = .016445;	/// Sherpa_CT10_ZqqZll_SHv21_improved
m_crossSectionMap[361095] = .006777;	/// Sherpa_CT10_WqqZvv_SHv21_improved
m_crossSectionMap[361094] = .0034234;	/// Sherpa_CT10_WqqZll_SHv21_improved
m_crossSectionMap[361093] = .011494;	/// Sherpa_CT10_WlvZqq_SHv21_improved
m_crossSectionMap[361092] = .024857;	/// Sherpa_CT10_WpqqWmlv_SHv21_improved
m_crossSectionMap[361091] = .024885;	/// Sherpa_CT10_WplvWmqq_SHv21_improved
m_crossSectionMap[361088] = .0034001;	/// Sherpa_CT10_lvvv
m_crossSectionMap[361086] = .01659;	/// Sherpa_CT10_ZqqZll
m_crossSectionMap[361084] = .0037583;	/// Sherpa_CT10_WqqZll
m_crossSectionMap[361068] = .014;	/// Sherpa_CT10_llvv
m_crossSectionMap[361067] = .0050248;	/// Sherpa_CT10_lllvOFPlus
m_crossSectionMap[361066] = .0025618;	/// Sherpa_CT10_lllvSFPlus
m_crossSectionMap[361065] = .0036254;	/// Sherpa_CT10_lllvOFMinus
m_crossSectionMap[361064] = .0018442;	/// Sherpa_CT10_lllvSFMinus
m_crossSectionMap[361063] = .012849;	/// Sherpa_CT10_llll
m_crossSectionMap[301534] = .024944;	/// Pythia8EvtGen_A14NNPDF23LO_Wprime_munu_Flat_13TeV
m_crossSectionMap[301533] = .02496;	/// Pythia8EvtGen_A14NNPDF23LO_Wprime_enu_Flat_13TeV
m_crossSectionMap[301138] = .000000000061624;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_5000M
m_crossSectionMap[301137] = .00000000014279;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_4500M5000
m_crossSectionMap[301136] = .0000000004721;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_4000M4500
m_crossSectionMap[301135] = .0000000015975;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_3500M4000
m_crossSectionMap[301134] = .0000000056692;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_3000M3500
m_crossSectionMap[301133] = .0000000072592;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2750M3000
m_crossSectionMap[301132] = .000000014549;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2500M2750
m_crossSectionMap[301131] = .000000030167;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2250M2500
m_crossSectionMap[301130] = .000000065311;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_2000M2250
m_crossSectionMap[301129] = .00000014945;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1750M2000
m_crossSectionMap[301128] = .00000036812;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1500M1750
m_crossSectionMap[301127] = .0000010029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1250M1500
m_crossSectionMap[301126] = .0000031594;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_1000M1250
m_crossSectionMap[301125] = .0000082865;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_800M1000
m_crossSectionMap[301124] = .00003098;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_600M800
m_crossSectionMap[301123] = .00017541;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_400M600
m_crossSectionMap[301122] = .0010832;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_250M400
m_crossSectionMap[301121] = .0032853;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_180M250
m_crossSectionMap[301120] = .022198;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminmunu_120M180
m_crossSectionMap[301118] = .00000000015346;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_5000M
m_crossSectionMap[301117] = .00000000040127;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_4500M5000
m_crossSectionMap[301116] = .0000000014307;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_4000M4500
m_crossSectionMap[301115] = .0000000050968;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_3500M4000
m_crossSectionMap[301114] = .00000001845;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_3000M3500
m_crossSectionMap[301113] = .000000023476;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2750M3000
m_crossSectionMap[301112] = .000000046259;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2500M2750
m_crossSectionMap[301111] = .000000093349;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2250M2500
m_crossSectionMap[301110] = .00000019463;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_2000M2250
m_crossSectionMap[301109] = .00000042457;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1750M2000
m_crossSectionMap[301108] = .00000098628;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1500M1750
m_crossSectionMap[301107] = .0000025071;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1250M1500
m_crossSectionMap[301106] = .0000072895;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_1000M1250
m_crossSectionMap[301105] = .000017668;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_800M1000
m_crossSectionMap[301104] = .000060793;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_600M800
m_crossSectionMap[301103] = .00031235;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_400M600
m_crossSectionMap[301102] = .0017543;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_250M400
m_crossSectionMap[301101] = .0050029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_180M250
m_crossSectionMap[301100] = .032053;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusmunu_120M180
m_crossSectionMap[301098] = .000000000061624;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_5000M
m_crossSectionMap[301097] = .00000000014279;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_4500M5000
m_crossSectionMap[301096] = .0000000004721;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_4000M4500
m_crossSectionMap[301095] = .0000000015975;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_3500M4000
m_crossSectionMap[301094] = .0000000056692;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_3000M3500
m_crossSectionMap[301093] = .0000000072592;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2750M3000
m_crossSectionMap[301092] = .000000014549;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2500M2750
m_crossSectionMap[301091] = .000000030167;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2250M2500
m_crossSectionMap[301090] = .000000065311;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_2000M2250
m_crossSectionMap[301089] = .00000014945;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1750M2000
m_crossSectionMap[301088] = .00000036812;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1500M1750
m_crossSectionMap[301087] = .0000010029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1250M1500
m_crossSectionMap[301086] = .0000031594;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_1000M1250
m_crossSectionMap[301085] = .0000082865;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_800M1000
m_crossSectionMap[301084] = .00003098;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_600M800
m_crossSectionMap[301083] = .00017541;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_400M600
m_crossSectionMap[301082] = .0010832;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_250M400
m_crossSectionMap[301081] = .0032852;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_180M250
m_crossSectionMap[301080] = .022198;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminenu_120M180
m_crossSectionMap[301078] = .00000000015346;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_5000M
m_crossSectionMap[301077] = .00000000040127;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_4500M5000
m_crossSectionMap[301076] = .0000000014307;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_4000M4500
m_crossSectionMap[301075] = .0000000050968;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_3500M4000
m_crossSectionMap[301074] = .00000001845;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_3000M3500
m_crossSectionMap[301073] = .000000023476;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2750M3000
m_crossSectionMap[301072] = .000000046259;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2500M2750
m_crossSectionMap[301071] = .000000093349;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2250M2500
m_crossSectionMap[301070] = .00000019463;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_2000M2250
m_crossSectionMap[301069] = .0000004245;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1750M2000
m_crossSectionMap[301068] = .00000098628;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1500M1750
m_crossSectionMap[301067] = .0000025071;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1250M1500
m_crossSectionMap[301066] = .0000072895;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_1000M1250
m_crossSectionMap[301065] = .000017668;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_800M1000
m_crossSectionMap[301064] = .000060793;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_600M800
m_crossSectionMap[301063] = .00031235;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_400M600
m_crossSectionMap[301062] = .0017543;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_250M400
m_crossSectionMap[301061] = .0050029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_180M250
m_crossSectionMap[301060] = .032053;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu_120M180
m_crossSectionMap[301038] = .00000000012649;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_5000M
m_crossSectionMap[301037] = .00000000028071;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_4500M5000
m_crossSectionMap[301036] = .00000000089764;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_4000M4500
m_crossSectionMap[301035] = .0000000029342;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_3500M4000
m_crossSectionMap[301034] = .000000010029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_3000M3500
m_crossSectionMap[301033] = .000000012487;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2750M3000
m_crossSectionMap[301032] = .000000024452;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2500M2750
m_crossSectionMap[301031] = .0000000494;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2250M2500
m_crossSectionMap[301030] = .00000010387;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_2000M2250
m_crossSectionMap[301029] = .00000022991;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1750M2000
m_crossSectionMap[301028] = .00000054521;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1500M1750
m_crossSectionMap[301027] = .0000014219;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1250M1500
m_crossSectionMap[301026] = .0000042582;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_1000M1250
m_crossSectionMap[301025] = .000010607;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_800M1000
m_crossSectionMap[301024] = .000037399;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_600M800
m_crossSectionMap[301023] = .0001955;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_400M600
m_crossSectionMap[301022] = .001082;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_250M400
m_crossSectionMap[301021] = .0029212;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_180M250
m_crossSectionMap[301020] = .017478;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYmumu_120M180
m_crossSectionMap[301018] = .00000000012649;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_5000M
m_crossSectionMap[301017] = .00000000028071;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_4500M5000
m_crossSectionMap[301016] = .00000000089764;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_4000M4500
m_crossSectionMap[301015] = .0000000029342;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_3500M4000
m_crossSectionMap[301014] = .000000010029;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_3000M3500
m_crossSectionMap[301013] = .000000012487;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2750M3000
m_crossSectionMap[301012] = .000000024452;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2500M2750
m_crossSectionMap[301011] = .0000000494;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2250M2500
m_crossSectionMap[301010] = .00000010387;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_2000M2250
m_crossSectionMap[301009] = .00000022992;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1750M2000
m_crossSectionMap[301008] = .00000054521;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1500M1750
m_crossSectionMap[301007] = .0000014219;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1250M1500
m_crossSectionMap[301006] = .0000042582;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_1000M1250
m_crossSectionMap[301005] = .000010607;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_800M1000
m_crossSectionMap[301004] = .000037401;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_600M800
m_crossSectionMap[301003] = .0001955;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_400M600
m_crossSectionMap[301002] = .001082;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_250M400
m_crossSectionMap[301001] = .0029212;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_180M250
m_crossSectionMap[301000] = .017478;	/// PowhegPythia8EvtGen_AZNLOCTEQ6L1_DYee_120M180
  
  
}

DataSetInfo::~DataSetInfo(){}

double DataSetInfo::getnEventsPerSample(unsigned int datasetID){
        
  try{
    return m_nEventsMap[datasetID];
  }
  catch(const std::out_of_range& oor){
    std::cout<<"[ERROR]\t DataSetInfo::getnEventsPerSample: Failed to find nEvents for dataset"
    " ID " << datasetID;
  }
  return 0.0;
        
}

double DataSetInfo::getCrossSection(unsigned int datasetID){
        
  try{
    return m_crossSectionMap[datasetID]; /// WARNING do we need to scale xsec?
  }
  catch(const std::out_of_range& oor){
    std::cout<<"[ERROR]\t DataSetInfo::getCrossSection: Failed to find crossSection for dataset"
    " ID " << datasetID;
  }
  return 0.0;
        
}
