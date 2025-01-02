# Solvers return equivalent solutions with POLYMOD data

    Code
      epi_outcome_iterative
    Output
        demo_grp       susc_grp susceptibility group_size p_infected n_infected
      1   [0,20)      immunised           0.10    4933097 0.02979109   146962.3
      2  [20,40)      immunised           0.10    5508767 0.02279411   125567.4
      3      40+      immunised           0.10    9653720 0.01603771   154823.5
      4   [0,20) part-immunised           0.55    4933097 0.15324273   755961.2
      5  [20,40) part-immunised           0.55    5508767 0.11910648   656129.9
      6      40+ part-immunised           0.55    9653720 0.08508340   821371.3
      7   [0,20)    susceptible           1.00    4933097 0.26098610  1287469.7
      8  [20,40)    susceptible           1.00    5508767 0.20592640  1134400.6
      9      40+    susceptible           1.00    9653720 0.14928409  1441146.7

---

    Code
      epi_outcome_newton
    Output
        demo_grp       susc_grp susceptibility group_size p_infected n_infected
      1   [0,20)      immunised           0.10    4933097 0.02979101   146961.9
      2  [20,40)      immunised           0.10    5508767 0.02279407   125567.2
      3      40+      immunised           0.10    9653720 0.01603768   154823.3
      4   [0,20) part-immunised           0.55    4933097 0.15324236   755959.4
      5  [20,40) part-immunised           0.55    5508767 0.11910629   656128.8
      6      40+ part-immunised           0.55    9653720 0.08508329   821370.2
      7   [0,20)    susceptible           1.00    4933097 0.26098551  1287466.8
      8  [20,40)    susceptible           1.00    5508767 0.20592609  1134398.9
      9      40+    susceptible           1.00    9653720 0.14928389  1441144.8

# Newton solver is correct in complex case

    Code
      epi_outcome
    Output
          demo_grp   susc_grp susceptibility group_size p_infected n_infected
      1 demo_grp_1 susc_grp_1              1   10831795  0.4656377  5043691.8
      2 demo_grp_2 susc_grp_1              1   11612456  0.4636757  5384413.8
      3 demo_grp_3 susc_grp_1              1   13511496  0.3987979  5388355.9
      4 demo_grp_4 susc_grp_1              1   11499398  0.3364829  3869350.8
      5 demo_grp_5 susc_grp_1              1    8167102  0.2469828  2017133.5
      6 demo_grp_6 susc_grp_1              1    4587765  0.1485668   681589.6

# Iterative solver is correct in complex case

    Code
      epi_outcome
    Output
          demo_grp   susc_grp susceptibility group_size p_infected n_infected
      1 demo_grp_1 susc_grp_1              1   10831795  0.4656377  5043692.4
      2 demo_grp_2 susc_grp_1              1   11612456  0.4636758  5384414.6
      3 demo_grp_3 susc_grp_1              1   13511496  0.3987979  5388355.7
      4 demo_grp_4 susc_grp_1              1   11499398  0.3364829  3869350.4
      5 demo_grp_5 susc_grp_1              1    8167102  0.2469827  2017133.2
      6 demo_grp_6 susc_grp_1              1    4587765  0.1485668   681589.4

