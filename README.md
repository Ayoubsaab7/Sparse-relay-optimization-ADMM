# About this repository
This repo contains files pertaining to my Master's thesis project, which was a communication system model under the umbrella of a Cloud Radio Access Network (C-RAN). The aforementioned files include:

1) The thesis
2) The IEEE-VTC-2018 (Chicago) Conference Paper: https://ieeexplore.ieee.org/document/8690755
3) Code

# Note to the reader
To understand the provided `C++` code, one must first go through the equations in the conference paper and play around with the math to get a sense of what's going on. However, the paper's abstract can be difficult to grasp if you are unfamiliar with the domain-specific terms. 

In what follows, I will do my best to 1) explain the motivation behind this project and 2) describe the topology of the communication system in layman terms.

# 1) Project Motivation
The performance of wireless networks has dramatically improved in the past 30 years, as transmission rates have risen by a thousand-fold from the first generation (1G) to the fourth generation (4G). Despite this feat, current networks are facing the challenge of meeting the supply for the projected data volumes. In fact, the proliferation of wireless devices with advanced capabilities has resulted in an unprecedented increase in cell density and -- more importantly -- throughput per user. Additionally, the emergence of new communication paradigms (e.g. vehicle-to-vehicle, device-to-device) has imposed significant challenges on the design and implementation of the next generation cellular networks (5G) as they aim to reduce:
 - end-to-end latency 
 - power consumption 
 - processing complexity 
 - cost 
5G is required to enable novel (and potentially conflicting) applications, making disruptive technologies and architectures indispensable for its genesis. 

To address the aforementioned data challenge and increase spectrum and energy efficiency, small-cell networks that complement traditional radio access networks (RAN) have been recognized as a potential solution. Specifically, centralized cloud-processing that is coupled with a dense deployment of low-complexity access points has been the discussion of many works. In this new architecture (termed C-RAN), the traditional base station (BS) functionalities are **apportioned** between the centralized datacenter and the access points. The former handles baseband signal processing and is thus termed the baseband unit (BBU) pool. The latter are responsible for data transmission/reception and analog-to-digital conversion and are termed remote radio heads (RRH) due to their distributed nature. The BBU and RRHs exchange channel state information (CSI) and user-equipment (UE) traffic data via low-latency and high-bandwidth optical transport links.

# 2) Topology and brief description
The communication system in question consists of `K` single-antenna user-equipments (e.g. mobile devices), each trying to communicate with another user-equipment. Accordingly, there are `K` pairs of source-destination single-antenna users. 

The wireless channels between these pairs are distorted but there is something that can be done to mitigate this issue. Namely, we place many multiple-antenna (MIMO) relay stations between them that amplify-and-forward (AF) the signals they receive from the sources. The received signal at each relay can be mathematically represented as a vector of size `N[l]`, where `N[l]` represents the number of antennas at relay `l`. 

The amplification step to be performed by each relay is simply a multiplication of the received vector by an appropriate matrix. This calls for the design of a matrix **for each relay** while satisfying some imposed constraints, such as transmit power constraint per relay.

This leads us to the following: 
 - Q. how should we amplify the signal? Alternatively, how do we design an appropriate matrix?   
 - A. It all depends on the chosen performance criterion. We usually either try to minimize a cost crieterion or maximize a benefit crietrion.

In this project, we choose to minimize a quantity called the total interference leakage. Therefore, the design of the amplification matrices is guided towards achieving this goal; again subject to power constraints for each relay and -- additionally -- constraints which preserve the desired signals at the destinations. In addition to the design of each matrix, we would like to save energy in this network by deactivating the ineffective relays. This calls for the modification of the chosen performance criterion to incorporate this requirement, which is termed as **regularization** in mathematical lingo.

# `C++` simulation
To experiment with the code: 
1) Download the "Sparse-relay-optimization-ADMM" folder cointaing all the code.
2) In `parameters.h`, set the variable `monteCarlo` to `1`.
3) Uncomment the line `//mimoObject.simulateTxRx(solution_vector);` in `main.cpp`.
4) On the terminal, compile `main.cpp` via the command: `g++ main.cpp -o experiment`
5) On the terminal, run the executable via the command: `./experiment`


# Contact: suggestions, comments, questions are welcome.
ayoub.saab@mail.mcgill.ca
