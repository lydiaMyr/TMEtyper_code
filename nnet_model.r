load("cluster_assign_new_tumor.rda")
load("TME_sub_exp.rda")
load("panCancer_cluster_data.rda")
load("immune_path_sig_top20_final.rda")

library(GSVA)
exp_mt = as.matrix(log2(exp_data+1))

ssgsea_param <- ssgseaParam(exprData = exp_mt,
                           geneSets = final_path_list,
                           alpha = 0.25,      # 排名权重指数
                           normalize = TRUE,       # 调节参数
                           minSize = 10,      # 基因集最小大小
                           maxSize = 500)     # 基因集最大大小

# 运行SSGSEA
ssgsea_results <- gsva(ssgsea_param, verbose = TRUE)



ft = row.names(pahtway_nes_all)
final_feature = final_path_list[ft]

train_data = t(pahtway_nes_all[ft,])
sample_label = cluster_assign_new[intersect(names(cluster_assign_new),colnames(pahtway_nes_all))]
test_data = t(ssgsea_results[ft,])

ft_new <- as.vector(unlist(lapply(ft, function(x) gsub("[ ]", "_", x))))
ft_new <- as.vector(unlist(lapply(ft_new, function(x) gsub("[/]", "_", x))))
ft_new <- as.vector(unlist(lapply(ft_new, function(x) gsub("[-]", "_", x))))
colnames(train_data) = ft_new
colnames(test_data) = ft_new


final_path_list_filter=final_path_list[ft]
names(final_path_list_filter)=ft_new
final_feature=final_path_list_filter
# 加载必要的包
library(nnet)    # 神经网络
library(caret)   # 数据预处理和模型评估
sample_label[which(sample_label=="TME1")]="MS"
sample_label[which(sample_label=="TME2")]="APH"
sample_label[which(sample_label=="TME3")]="LRH"
sample_label[which(sample_label=="TME4")]="KSH"
sample_label[which(sample_label=="TME5")]="MCH"
sample_label[which(sample_label=="TME6")]="FSH"
sample_label[which(sample_label=="TME7")]="MSH"
sample_label <- as.factor(sample_label)
df=data.frame(Sample=names(sample_label),Subtype=sample_label)
# 合并训练数据和标签
train_df <- data.frame(train_data[names(sample_label),], Class = sample_label)
train_data=train_data[names(sample_label),]
## 数据预处理
# 标准化数据(中心化和缩放)
preProc <- preProcess(train_data, method = c("center", "scale"))
train_data_std <- predict(preProc, train_data)
test_data_std <- predict(preProc, test_data)

## 训练神经网络模型
# 设置随机种子保证可重复性
set.seed(123)

#神经网络模型
model <- nnet(Class ~ ., 
              data = data.frame(train_data_std, Class = sample_label),
              size = 64,           # 增加隐藏层神经元数量（匹配CNN的滤波器数量）
              maxit = 100,         # 减少迭代次数（配合早停机制）
              MaxNWts = 50000,     # 增加最大权重数量（适应更大的网络结构）
              decay = 0.001,       # 减小权重衰减（匹配Adam优化器的学习率）
              abstol = 1.0e-6,     # 添加绝对容忍度（早停条件）
              reltol = 1.0e-8,     # 添加相对容忍度（早停条件）
              trace = TRUE)        # 显示训练过程
save(model,file="nnet_model.rda")
## 模型评估(训练集)
train_pred <- predict(model, train_data_std, type = "class")
confusionMatrix(as.factor(train_pred), sample_label)

## 预测新数据
test_pred1 <- predict(model, test_data_std, type = "class")
test_prob1 <- predict(model, test_data_std, type = "raw")

# 输出预测结果
print("测试数据的预测类别:")
print(test_pred1)

# 输出预测概率(7个类别的概率)
print("测试数据的预测概率:")
print(test_prob1)

