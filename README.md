# GRN-inference-combining-scJoint-and-Pando
本软件的运行环境为R 4.1.2和python 3.7.10。完整程序包括demo数据、源代码、公共库以及产生的结果文件已上传  
demo_data文件夹：存放的是人类外周血RNA和ATAC测试demo数据，由于大小超过100M.，因此将其分卷压缩，使用时需要包含分卷在内一起解压；  
data文件夹：存放的是数据预处理得到的结果；  
models文件夹：scJoint运行后产生的训练模型；  
output文件夹：scJoint运行后得到的输出结果；  
scJoint_result文件夹：将输出结果进行整合为pando输入做准备；  
pando_result文件夹：最终结果；  
util文件夹：必需的依赖工具、函数或公共代码。  
总程序代码主要包括六个步骤的集成，分别是：  
第一步：data_to_h5.R 数据预处理，将表达矩阵转成h5文件；  
第二步：h5_to_npz.py 数据预处理，将h5文件转成npz文件；  
第三步：config.py 配置文件；  
第四步：scjoint.py 整合RNA-seq数据和ATAC-seq数据程序主函数；  
第五步：umap_embedding.R scjoint处理后得到的结果处理；  
第六步：pando.R 推断基因调控网络。  
以及还包括一些必需的依赖工具、函数或公共代码存放在util文件夹下。  
