# GRN-inference-combining-scJoint-and-Pando
## 软件简介
利用**注释好的RNA-seq数据**，来实现将细胞类型的标记信息合并到一个半监督范式中来训练**无标记的ATAC-seq数据**，最后得到ATAC-seq数据细胞类型的预测身份以及联合降维信息。再利用上述得到的联合降维矩阵信息将RNA和ATAC数据的细胞进行配对，并将两个数据集的细胞barcode重命名，使其两个数据集的细胞barcode名字一致，组合成一个包含双组学的pseudo multi对象，然后利用pando对其进行基因调控网络推断。
## 软件支持及组成
本软件的运行环境为R 4.1.2和python 3.7.10，此外Pando包最新1.0.5版本会报错，建议使用1.0.3版本。完整程序包括demo数据、源代码、公共库以及产生的结果文件已上传  
__demo_data文件夹__：存放的是人类外周血RNA和ATAC测试demo数据，由于大小超过100M.，因此将其分卷压缩，使用时需要包含分卷在内一起解压；  
__data文件夹__：存放的是数据预处理得到的结果；  
__models文件夹__：scJoint运行后产生的训练模型；  
__output文件夹__：scJoint运行后得到的输出结果；  
__scJoint_result文件夹__：将输出结果进行整合为pando输入做准备；  
__pando_result文件夹__：最终结果；  
__util文件夹__：scJoint运行必需的依赖工具、函数或公共代码。  
总程序代码主要包括六个步骤的集成，分别是：  
__第一步__：data_to_h5.R 数据预处理，将表达矩阵转成h5文件；  
__第二步__：h5_to_npz.py 数据预处理，将h5文件转成npz文件；  
__第三步__：config.py 配置文件；  
__第四步__：scjoint.py 整合RNA-seq数据和ATAC-seq数据程序主函数；  
__第五步__：umap_embedding.R scjoint处理后得到的结果处理；  
__第六步__：pando.R 推断基因调控网络。  
## 主要结果
### scJoint整合UMAP图
![RNA_ATAC](https://github.com/tyc-1998/GRN-inference-combining-scJoint-and-Pando/assets/58178820/12924155-8e5c-4423-9964-3bf96412c7b8)
### pando推断GRN结果
![GRN_plot](https://github.com/tyc-1998/GRN-inference-combining-scJoint-and-Pando/assets/58178820/363a25bf-4e91-432f-a320-0b5ee2f2bb1b)

