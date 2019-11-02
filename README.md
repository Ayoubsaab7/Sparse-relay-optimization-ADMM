# About this repository
This repo contains files pertaining to my Master's thesis project, which was a communication system model under the umbrella of a Cloud Radio Access Network (C-RAN). The aforementioned files include:


1) The thesis
2) The IEEE-VTC-2018 (Chicago) Conference Paper: https://ieeexplore.ieee.org/document/8690755
3) C++ code

# Note to the reader
To understand the provided code, one must go through the conference paper and play around with the math. However, the paper's abstract can be difficult to grasp for people who don't know domain-specific terms. In what follows, I will do my best to explain the purpose of the project and the topology of the communication system in layman terms.

# Project Description
The communication system consists of 'K' single-antenna user-equipments (e.g. mobile phones), each trying to communicate with another user-equipment. Accordingly, there are 'K' pairs of source-destination single-antenna user-equipments. The wireless channel between these pairs is distorted; but there is something that can be done to mitigate this issue: we place many multiple-antenna (MIMO) relay stations that have as function to amplify and forward the signals they receive. For this reason, the relay stations are called amplify-and-forward (AF) relays to distinguish them from other types such as decode-and-forward (DF). 

The received signal at each relay can be mathematically represented as a vector of size 'NL', where 'NL' stands for the number of antennas at relay 'L'. The amplification step that each relay will perform is simply a multiplication of the received vector by an appropriate matrix. This leads us to the following question: how should we amplify the signal? Said differently, how do we pick a 'good' matrix? The design of an amplification matrix must be done for each relay and must satisfy some imposed constraints (e.g. power constraint). The answer to this question depends upon the chosen criterion that we usually either try to minimize or maximize. In this project, we have decided to minimize a quantity called the total interference leakage. Therefore, the design of the matrices is oriented towards achieving the aforementioned objective, again subject to power constraints and additionally, distortionless constraints (which preserve the desired signals at the destinations).

In addition to the design of each matrix, we would like to save some energy in this network configuration by turning off the ineffective relays. This can be mathematically done through a process called regularization.

# C++ simulation


# Contact: suggestions, comments, questions are welcome.
ayoub.saab@mail.mcgill.ca
