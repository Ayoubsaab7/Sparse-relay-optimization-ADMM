# About this repository
This repo contains files pertaining to my Master's thesis project, which was a communication system model under the umbrella of a Cloud Radio Access Network (C-RAN). The aforementioned files include:

1) The thesis
2) The IEEE-VTC-2018 (Chicago) Conference Paper: https://ieeexplore.ieee.org/document/8690755
3) C++ code

# Note to the reader
To understand the provided C++ code, one must first go through the equations in the conference paper and play around with the math to get sense of what's going on. However, the paper's abstract can be difficult to grasp for people who don't know domain-specific terms. In what follows, I will do my best to explain the motivation behind this project and describe the topology of the communication system in layman terms.

# Project Motivation
The performance of wireless networks has dramatically improved in the past 30 years, as transmission rates have risen by a thousand-fold from the first generation (1G) to the fourth generation (4G). Despite this feat, current networks are facing a great challenge as they will not be able to supply the projected data volumes. In fact, the proliferation of wireless devices with advanced capabilities has resulted in an unprecedented increase in both cell density and, more importantly, throughput per user. Additionally, the emergence of new communication paradigms (e.g. vehicle-to-vehicle, device-to-device) has imposed significant challenges on the design and implementation of the next generation cellular networks (5G) as they aim to reduce end-to-end latency, power consumption, processing complexity and cost. 5G is required to enable novel (and potentially conflicting) applications, making disruptive technologies and architectures indispensable for its genesis.

To address the aforementioned data challenge and increase spectrum and energy efficiency, small cell networks complementing traditional radio access networks (RAN) have been recognized as a potential solution. Specifically, centralized cloud-processing coupled with a dense deployment of low-complexity access points has been the discussion of many works. In this new architecture (termed C-RAN), the traditional base station (BS) functionalities are apportioned between the centralized datacenter and the access points. The former handles baseband signal processing functions and is thus termed the baseband unit (BBU) pool. The latter are responsible for data transmission/reception and analog-to-digital conversion and are termed remote radio heads (RRH) due to their distributed nature. The BBU and RRHs exchange channel state information (CSI) and user-equipment (UE) traffic data via low-latency and high-bandwidth optical transport links.

# Topology and description
The communication system consists of 'K' single-antenna user-equipments (e.g. mobile phones), each trying to communicate with another user-equipment. Accordingly, there are 'K' pairs of source-destination single-antenna user-equipments. The wireless channel between these pairs is distorted; but there is something that can be done to mitigate this issue: we place many multiple-antenna (MIMO) relay stations that have as function to amplify and forward the signals they receive. For this reason, the relay stations are called amplify-and-forward (AF) to distinguish them from other types such as decode-and-forward (DF). 

The received signal at each relay can be mathematically represented as a vector of size 'NL', where 'NL' stands for the number of antennas at relay 'L'. The amplification step that each relay will perform is simply a multiplication of the received vector by an appropriate matrix. This leads us to the following question: how should we amplify the signal? Said differently, how do we pick a 'good' matrix? The design of an amplification matrix must be done for each relay and must satisfy some imposed constraints (e.g. power constraint). The answer to this question depends upon the chosen criterion that we usually either try to minimize or maximize. In this project, we have decided to minimize a quantity called the total interference leakage. Therefore, the design of the matrices is oriented towards achieving the aforementioned objective, again subject to power constraints and additionally, distortionless constraints (which preserve the desired signals at the destinations).

In addition to the design of each matrix, we would like to save some energy in this network configuration by turning off the ineffective relays. This can be mathematically done through a process called regularization.

# C++ simulation
To run the code, download the "Sparse-relay-optimization-ADMM" folder cointaing "main.cpp" and "Eigen".

# Contact: suggestions, comments, questions are welcome.
ayoub.saab@mail.mcgill.ca
