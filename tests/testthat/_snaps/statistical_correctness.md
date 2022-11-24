# Solvers return equivalent solutions with POLYMOD data

    Code
      epi_outcome_iterative
    Output
        demo_grp       susc_grp susceptibility p_infected
      1   [0,20)      immunised           0.10 0.02979109
      2  [20,40)      immunised           0.10 0.02279411
      3      40+      immunised           0.10 0.01603771
      4   [0,20) part-immunised           0.55 0.15324273
      5  [20,40) part-immunised           0.55 0.11910648
      6      40+ part-immunised           0.55 0.08508340
      7   [0,20)    susceptible           1.00 0.26098610
      8  [20,40)    susceptible           1.00 0.20592640
      9      40+    susceptible           1.00 0.14928409

---

    Code
      epi_outcome_newton
    Output
        demo_grp       susc_grp susceptibility p_infected
      1   [0,20)      immunised           0.10 0.02979101
      2  [20,40)      immunised           0.10 0.02279407
      3      40+      immunised           0.10 0.01603768
      4   [0,20) part-immunised           0.55 0.15324236
      5  [20,40) part-immunised           0.55 0.11910629
      6      40+ part-immunised           0.55 0.08508329
      7   [0,20)    susceptible           1.00 0.26098551
      8  [20,40)    susceptible           1.00 0.20592609
      9      40+    susceptible           1.00 0.14928389

# Newton solver is correct in complex case

    Code
      epi_outcome
    Output
          demo_grp   susc_grp susceptibility p_infected
      1 demo_grp_1 susc_grp_1              1  0.4656377
      2 demo_grp_2 susc_grp_1              1  0.4636757
      3 demo_grp_3 susc_grp_1              1  0.3987979
      4 demo_grp_4 susc_grp_1              1  0.3364829
      5 demo_grp_5 susc_grp_1              1  0.2469828
      6 demo_grp_6 susc_grp_1              1  0.1485668

# Iterative solver is correct in complex case

    Code
      epi_outcome
    Output
          demo_grp   susc_grp susceptibility p_infected
      1 demo_grp_1 susc_grp_1              1  0.4656377
      2 demo_grp_2 susc_grp_1              1  0.4636758
      3 demo_grp_3 susc_grp_1              1  0.3987979
      4 demo_grp_4 susc_grp_1              1  0.3364829
      5 demo_grp_5 susc_grp_1              1  0.2469827
      6 demo_grp_6 susc_grp_1              1  0.1485668

