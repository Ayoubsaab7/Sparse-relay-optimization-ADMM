# About this repository
This repo contains files pertaining to my Master's thesis project, which was a communication system model under the umbrella of a Cloud Radio Access Network (C-RAN). The aforementioned files include:

1) The thesis
2) The IEEE-VTC-2018 (Chicago) Conference Paper: https://ieeexplore.ieee.org/document/8690755
3) C++ code

# Note to the reader
To understand the provided C++ code, one must first go through the equations in the conference paper and play around with the math to get a sense of what's going on. However, the paper's abstract can be difficult to grasp for people who don't know domain-specific terms. In what follows, I will do my best to 1) explain the motivation behind this project and 2) describe the topology of the communication system in layman terms.

# 1) Project Motivation
The performance of wireless networks has dramatically improved in the past 30 years, as transmission rates have risen by a thousand-fold from the first generation (1G) to the fourth generation (4G). Despite this feat, current networks are facing a great challenge as they will not be able to supply the projected data volumes. In fact, the proliferation of wireless devices with advanced capabilities has resulted in an unprecedented increase in both cell density and, more importantly, throughput per user. Additionally, the emergence of new communication paradigms (e.g. vehicle-to-vehicle, device-to-device) has imposed significant challenges on the design and implementation of the next generation cellular networks (5G) as they aim to reduce end-to-end latency, power consumption, processing complexity and cost. 5G is required to enable novel (and potentially conflicting) applications, making disruptive technologies and architectures indispensable for its genesis.

To address the aforementioned data challenge and increase spectrum and energy efficiency, small cell networks complementing traditional radio access networks (RAN) have been recognized as a potential solution. Specifically, centralized cloud-processing coupled with a dense deployment of low-complexity access points has been the discussion of many works. In this new architecture (termed C-RAN), the traditional base station (BS) functionalities are apportioned between the centralized datacenter and the access points. The former handles baseband signal processing functions and is thus termed the baseband unit (BBU) pool. The latter are responsible for data transmission/reception and analog-to-digital conversion and are termed remote radio heads (RRH) due to their distributed nature. The BBU and RRHs exchange channel state information (CSI) and user-equipment (UE) traffic data via low-latency and high-bandwidth optical transport links.

# 2) Topology and brief description
The communication system in question consists of 'K' single-antenna user-equipments (e.g. mobile devices), each trying to communicate with another user-equipment. Accordingly, there are 'K' pairs of source-destination single-antenna user-equipments. The wireless channels between these pairs are distorted but there is something that can be done to mitigate this issue. Namely, we place many multiple-antenna (MIMO) relay stations in between that amplify and forward the signals they receive from the sources. For this reason, the relay stations are called amplify-and-forward (AF) to distinguish them from other variations (such as decode-and-forward). 

The received signal at each relay can be mathematically represented as a vector of size 'NL', where 'NL' represents the number of antennas at relay 'L'. The amplification step to be performed by each relay is simply a multiplication of the received vector by an appropriate matrix. This leads us to the following question: how should we amplify the signal? Said differently, how do we design an appropriate matrix? This design should be done for each single relay and must satisfy some imposed constraints (e.g. transmitted power constraint).  

The answer to the previous questions depends on the chosen performance criterion that we usually either try to minimize (cost) or maximize (benefit). This calls for the employment of tools from mathematical optimization. In this project, we choose to minimize a quantity called the total interference leakage. Therefore, the design of the amplification matrices is guided towards achieving the aforementioned goal; again subject to power constraints for each relay and -- additionally -- distortionless constraints which preserve the desired signals at the destinations. In addition to the design of each matrix, we would like to save energy in this network configuration by deactivating the ineffective relays. This naturally calls for the modification of the chosen performance criterion to incorporate this extra requirement, which is termed as regularization in mathematical lingo.

# C++ simulation
To run the code, download the "Sparse-relay-optimization-ADMM" folder cointaing "main.cpp", "parameters.h", "functions.h" and "Eigen".

# Contact: suggestions, comments, questions are welcome.
ayoub.saab@mail.mcgill.ca
