clear 
Y=2000+[4 ...butyrate
    7  ...molecular epidemiology
    8  ...five human
    8  ...different
    8  ...epidemiology
    9  ...functional
    9  ...three novel
    10 ...congenital
    11 ...high prevalence
    12 ...extensive
    13 ...Erythrocytosis
    14 ...Intracranial
    17 ...Anthropometric
    17 ...twelve
    20 ...Clotting
    20 ...Fondaparinux
    20 ...Thromboprofilaxys
    20 ...Prognostic
    21 ...Antithrombotic
    21 ...Clinical differences
    21 ...The impact of risk
    21 ...Possible
    21 ...Pathophysiology 
    23 ...LUS Score
    23 ...PLASMIC
    26 ...Globalization
    26 ...Utility
    26 ...Beyond
    26 ...Clinical spectrum
    26 ...Papilledema
    ];
A=     [8  ... butyrate
    9  ...molecular epidemiology
    10 ...five human
    7  ...different
    11 ...epidemiology
    11 ...functional
    7  ...three novel 
    5  ...congenital 
    4  ...high prevalence
    8  ...extensive
    8  ...Erythrocytosis
    8  ...Intracranial 
    8  ...Anthropometric
    14 ...twelve
    7  ...Clotting
    13 ...Fondaparinux
    14 ...Thromboprofilaxys
    6  ...Prognostic
    14 ...Antithormbotic
    8 ...Clinical differences
    14 ...The impact of risk
    3 ...Possible
    9 ...Pathophysiology   
    16 ...LUS Score
    8 ...PLASMIC
    8 ...Globalization
    6 ...Utility
    1 ...Beyond
    7 ...Clinical spectrum
    14 ...Papilledema
    ];
C=     [81 ... butyrate
        32 ...molecular epidemiology
        13 ...five human
        2 ...different
        54 ...epidemiology
        31 ...functional
        35 ...three novel
        62 ...congenital
        20 ...high prevalence
        56 ...extensive
        1  ...Erythrocytosis
        49 ...Intracranial
        5  ...Anthropometric
        9  ...twelve
        61 ...Clotting
        25 ...Fondaparinux
        21 ...Thromboprofilaxys
        29 ...Prognostic
        25 ...Antithrombotic
        6 ...Clinical differences
        9 ...The impact of risk
        22 ...Possible
        10 ...Pathophysiology    
        0 ...LUS Score
        5 ...PLASMIC
        0 ...Globalization
        2 ...Utility
        1 ...Beyond
        1 ...Clinical spectrum
        0 ...Papilledema
    ];

R=bibliometrics(C,Y,A);
clear A Y C