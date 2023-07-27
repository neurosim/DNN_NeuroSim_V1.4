# DNN+NeuroSim V1.4

The DNN+NeuroSim framework was developed by [Prof. Shimeng Yu's group](https://shimeng.ece.gatech.edu/) (Georgia Institute of Technology). The model is made publicly available on a non-commercial basis. Copyright of the model is maintained by the developers, and the model is distributed under the terms of the [Creative Commons Attribution-NonCommercial 4.0 International Public License](http://creativecommons.org/licenses/by-nc/4.0/legalcode)

:star2: This is the released version 1.4 (August 1, 2023) for the tool, and this version has **_improved following inference engine estimation_**:
```
1. Support for technology scaling down to 1nm node in the hardware estimation framework (C++ code).

The following is a list of the supported nodes:

130nm
90nm
65nm
45nm
32nm
22nm
14nm
10nm
7nm
5nm   -   IRDS Projection (2021-2022)
3nm   -   IRDS Projection (2021-2022)
2nm   -   IRDS Projection (2021-2022)
1nm   -   IRDS Projection (2021-2022)

```

:point_right: :point_right: :point_right: **In "Param.cpp", to pick technology node:**
```
modify the "tech" parameter (line 159 in Param.cpp) to match the desired case
For example tech = 5 corresponds to a technode of 22nm and tech = 6 corresponds to 14nm.

For additional details about the device parameters used in NeuroSim, refer to section 7 of the V1.4 manual.
```

2. Add partial parallel mode in python wrapper (for single-level cells only) and C++ code for hardware estimation.

:point_right: :point_right: :point_right: **To enable partial parallel mode**
```
Specify the parameter "parallelRead" when running the python wrapper.
Partial parallel mode will be enabled in both the python wrapper and C++ code.

--parallelRead N

Where N is the desired number of rows activated in parallel and N <= sub-array size.
```

3. XY Bus as an alternative to H-tree interconnect.

:point_right: :point_right: :point_right: **In "Param.cpp", to switch interconnect mode:**
```
globalBusType = false;		// false: X-Y Bus      // true: H-Tree
```


**_three default examples for quick start_**:
```
1. VGG8 on cifar10 
   8-bit "WAGE" mode pretrained model is uploaded to './log/VGG8.pth'
3. DenseNet40 on cifar10 
   8-bit "WAGE" mode pretrained model is uploaded to './log/DenseNet40.pth'
5. ResNet18 on imagenet 
   "FP" mode pretrained model is loaded from 'https://download.pytorch.org/models/resnet18-5c106cde.pth'
```
:point_right: :point_right: :point_right: **To quickly start inference estimation of default models (skip training)**
```
python inference.py --dataset cifar10 --model VGG8 --mode WAGE --cellBit 1 --subArray 128 --parallelRead 64
python inference.py --dataset cifar10 --model DenseNet40 --mode WAGE
python inference.py --dataset imagenet --model ResNet18 --mode FP
```

<br/>

**_For estimation of on-chip training accelerators, please visit released V2.1 [DNN+NeuroSim V2.1](https://github.com/neurosim/DNN_NeuroSim_V2.1)_**

In Pytorch/Tensorflow wrapper, users are able to define **_network structures, precision of synaptic weight and neural activation_**. With the integrated NeuroSim which takes real traces from wrapper, the framework can support hierarchical organization from device level to circuit level, to chip level and to algorithm level, enabling **_instruction-accurate evaluation on both accuracy and hardware performance of inference_**.

Developers: [Junmo Lee](mailto:junmolee@gatech.edu) :two_men_holding_hands: [James Read](mailto:jread6@gatech.edu) :two_men_holding_hands: [Anni Lu](mailto:alu75@gatech.edu) :two_women_holding_hands: [Xiaochen Peng](mailto:xpeng76@gatech.edu) :two_women_holding_hands: [Shanshi Huang](mailto:shuang406@gatech.edu).

This research is supported by NSF CAREER award, NSF/SRC E2CDA program, PRISM center and CHIMES center, both part of the SRC/DARPA JUMP 2.0 program.

If you use the tool or adapt the tool in your work or publication, you are required to cite the following reference:

**_X. Peng, S. Huang, Y. Luo, X. Sun and S. Yu, ※[DNN+NeuroSim: An End-to-End Benchmarking Framework for Compute-in-Memory Accelerators with Versatile Device Technologies](https://ieeexplore-ieee-org.prx.library.gatech.edu/document/8993491), *§ IEEE International Electron Devices Meeting (IEDM)*, 2019._**

If you have logistic questions or comments on the model, please contact :man: [Prof. Shimeng Yu](mailto:shimeng.yu@ece.gatech.edu), and if you have technical questions or comments, please contact :man: [Junmo Lee](mailto:junmolee@gatech.edu) or :man: [James Read](mailto:jread6@gatech.edu) or :woman: [Anni Lu](mailto:alu75@gatech.edu).


## File lists
1. Manual: `Documents/DNN NeuroSim V1.4 Manual.pdf`
2. DNN_NeuroSim wrapped by Pytorch: 'Inference_pytorch'
3. NeuroSim under Pytorch Inference: 'Inference_pytorch/NeuroSIM'


## Installation steps (Linux + Anaconda/Miniconda)
We have included an Anaconda environment with this version to make package installation easier.

If you don't want to use the conda environment or don't have a CUDA enabled GPU, check the environment.yml file for the versions of all packages used.

This version supports the recently released PyTorch 2.0

We have currently tested the following CUDA drivers:


1. Download Anaconda/Miniconda: https://docs.conda.io/en/latest/miniconda.html
2. Follow install instructions: https://docs.conda.io/en/latest/miniconda.html#installing

3. Get the tool from GitHub
```
git clone https://github.com/neurosim/DNN_NeuroSim_V1.4.git
cd DNN_NeuroSim_V1.4
```

4. Create conda environemnt from provided environment file

```
conda create env --file environment.yml
```

5. Activate neurosim environment

```
conda activate neurosim
```

6. Train the network to get the model for inference (can be skipped by using pretrained default models)

7. Compile the NeuroSim codes
```
cd Inference_pytorch/NeuroSIM
make
```

8. Run Pytorch/Tensorflow wrapper with default settings (integrated with NeuroSim)

```
cd ..
python inference.py
```


For additional details on the usage of this tool, please refer to the manual.


## References related to this tool 
1. A. Lu, X. Peng, W. Li, H. Jiang, S. Yu, ※NeuroSim simulator for compute-in-memory hardware accelerator: validation and benchmark, *§ Frontiers in Artificial Intelligence, vol. 4, 659060, 2021.
2. X. Peng, S. Huang, Y. Luo, X. Sun and S. Yu, ※DNN+NeuroSim: An End-to-End Benchmarking Framework for Compute-in-Memory Accelerators with Versatile Device Technologies, *§ IEEE International Electron Devices Meeting (IEDM)*, 2019.
3. X. Peng, R. Liu, S. Yu, ※Optimizing weight mapping and data flow for convolutional neural networks on RRAM based processing-in-memory architecture, *§ IEEE International Symposium on Circuits and Systems (ISCAS)*, 2019.
4. P.-Y. Chen, S. Yu, ※Technological benchmark of analog synaptic devices for neuro-inspired architectures, *§ IEEE Design & Test*, 2019.
5. P.-Y. Chen, X. Peng, S. Yu, ※NeuroSim: A circuit-level macro model for benchmarking neuro-inspired architectures in online learning, *§ IEEE Trans. CAD*, 2018.
6. X. Sun, S. Yin, X. Peng, R. Liu, J.-S. Seo, S. Yu, ※XNOR-RRAM: A scalable and parallel resistive synaptic architecture for binary neural networks,*§ ACM/IEEE Design, Automation & Test in Europe Conference (DATE)*, 2018.
7. P.-Y. Chen, X. Peng, S. Yu, ※NeuroSim+: An integrated device-to-algorithm framework for benchmarking synaptic devices and array architectures, *§ IEEE International Electron Devices Meeting (IEDM)*, 2017.
8. P.-Y. Chen, S. Yu, ※Partition SRAM and RRAM based synaptic arrays for neuro-inspired computing,*§ IEEE International Symposium on Circuits and Systems (ISCAS)*, 2016.
9. P.-Y. Chen, D. Kadetotad, Z. Xu, A. Mohanty, B. Lin, J. Ye, S. Vrudhula, J.-S. Seo, Y. Cao, S. Yu, ※Technology-design co-optimization of resistive cross-point array for accelerating learning algorithms on chip,*§ IEEE Design, Automation & Test in Europe (DATE)*, 2015.
10. S. Wu, et al., ※Training and inference with integers in deep neural networks,*§ arXiv: 1802.04680*, 2018.
11. github.com/boluoweifenda/WAGE
12. github.com/stevenygd/WAGE.pytorch
13. github.com/aaron-xichen/pytorch-playground
