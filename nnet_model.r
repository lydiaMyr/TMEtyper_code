load("F:/workspace/TME_project_file/cluster_assign_new_tumor.rda")
load("F:/workspace/TME_project_file/文章修改/TME_sub_exp.rda")
load("F:/workspace/TME_project_file/panCancer_cluster_data0923.rda")
load("F:/workspace/TME_project_file/immune_path_sig_top20_final.rda")
#pathway feature final_path_list
#final_path_list
load("F:/workspace/免疫平衡评分工具/dataset/GSE53419_tpm_final.rda")
library(GSVA)

exp_mt = as.matrix(log2(GSE53419_tpm_final+1))
ssgsea_scores <- gsva(exp_mt, 
                      final_path_list,
                      method = "ssgsea",
                      kcdf = "Gaussian",  # 对于log2转换后的数据使用Gaussian
                      verbose = TRUE)


final_feature = final_path_list[ft]
save(final_feature,file="F:/workspace/TME_project_file/nnet_feature.rda")
# 结果是一个矩阵 (gene sets × samples)
ft = intersect(row.names(pahtway_nes_all),row.names(ssgsea_scores))[-208]
train_data = t(pahtway_nes_all[ft,])
sample_label = cluster_assign_new[colnames(pahtway_nes_all)]
test_data = t(ssgsea_scores[ft,])

ft_new <- as.vector(unlist(lapply(ft, function(x) gsub("[ ]", "_", x))))
ft_new <- as.vector(unlist(lapply(ft_new, function(x) gsub("[/]", "_", x))))
ft_new <- as.vector(unlist(lapply(ft_new, function(x) gsub("[-]", "_", x))))
colnames(train_data) = ft_new
colnames(test_data) = ft_new

# 加载必要的包
library(nnet)    # 神经网络
library(caret)   # 数据预处理和模型评估

## 数据准备
# 假设:
# train_data - 训练数据矩阵或数据框(行是样本，列是特征)
# sample_label - 训练数据的分类标签(7个类别)
# test_data - 测试数据(与train_data相同的特征)

# 确保数据格式正确
# 将标签转换为因子(7个水平)
sample_label <- as.factor(sample_label)

# 合并训练数据和标签
train_df <- data.frame(train_data, Class = sample_label)

## 数据预处理
# 标准化数据(中心化和缩放)
save(train_data,file="F:/workspace/TME_project_file/train_data.rda")
preProc <- preProcess(train_data, method = c("center", "scale"))
train_data_std <- predict(preProc, train_data)
test_data_std <- predict(preProc, test_data)

## 训练神经网络模型
# 设置随机种子保证可重复性
# set.seed(123)

# # 训练模型
# # size参数控制隐藏层神经元数量，可以根据数据复杂度调整
# # maxit控制最大迭代次数
# model <- nnet(Class ~ ., 
#               data = data.frame(train_data_std, Class = sample_label),
#               size = 10,       # 隐藏层神经元数量
#               maxit = 500,     # 最大迭代次数
#               MaxNWts = 10000,  # 最大权重数量
#               decay = 0.01,     # 权重衰减(正则化)
#               trace = TRUE)     # 显示训练过程
save(model,file="F:/workspace/TME_project_file/nnet_model.rda")
## 模型评估(训练集)
train_pred <- predict(model, train_data_std, type = "class")
confusionMatrix(as.factor(train_pred), sample_label)

## 预测新数据
test_pred <- predict(model, test_data_std, type = "class")
test_prob <- predict(model, test_data_std, type = "raw")

# 输出预测结果
print("测试数据的预测类别:")
print(test_pred)

# 输出预测概率(7个类别的概率)
print("测试数据的预测概率:")
print(test_prob)

