# About this repository
This repo contains files pertaining to my Master's thesis project, which was about a potential system model/scenario under the umbrella of a Cloud Radio Access Network (C-RAN). The aforementioned files
include:


0) The thesis
1) The IEEE-VTC-2018 (Chicago) Conference Paper: https://ieeexplore.ieee.org/document/8690755
2) C++ code
3) MATLAB code

# Note to the reader
The conference paper's abstract can be difficult to grasp for people who don't know domain-specific terms, the domain in question being Communication Systems. In what follows, I will do my best to explain the purpose of the project in layman terms. Let's begin.

# Project Description
We have a certain number (K) of mobile phones, every one of which paired with a single phone, with which it is trying to communicate. Unfortunately, the channel in between the K pairs is distorted. 
There is something that can be done to mitigate this issue: we place many multiple-antenna (MIMO) relays stations that have as function 
to first 1) amplify the signal they receive and second 2) forward it to the mobile stations that are expecting messages. For this reason
the relays are called amplify-and-forward (AF) relays.

Now mathematically, the received signal at relay 'l' is a vector of size 'N'. For instance, relay #1 has 2 antennas so the signal it 
receives is a vector of size 2. The amplification step that each relay will perform is simply a multiplication of the received vector by a matrix. 
This leads us to the following question: How should we amplify the signal? Said differently, how do we pick a 'good' matrix? 
This question must be answered for each relay. Finally, the design of each matrix must satisfy some imposed constraints (such as power constraints).

In addition to the design of each matrix, we would like to save some energy in this network configuration by turning off some of the 'useless' relays.
This process can be done using what is called regularization.

# C++ simulation

# MATLAB simulation

# Contact: suggestions, comments, questions are welcome.
ayoub.saab@mail.mcgill.ca
