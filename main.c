#define _USE_MATH_DEFINES

#include<stdio.h>
#include<inttypes.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<omp.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

// 伪随机数生成器

typedef struct {
    uint64_t state;
    uint64_t inc;
} pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng) {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// TODO 外置
#define HEAD "{\"header\":{\"layout\":10,\"icons\":[0,0,0,0,0],\"time\":\"2023-02-29T24:00:01.000Z\",\"gameVersion\":\"0.9.27.15466\",\"shortDesc\":\"New Blueprint\",\"desc\":\"\"},\"version\":1,\"cursorOffset\":{\"x\":0,\"y\":0},\"cursorTargetArea\":0,\"dragBoxSize\":{\"x\":1,\"y\":1},\"primaryAreaIdx\":0,\"areas\":[{\"index\":0,\"parentIndex\":-1,\"tropicAnchor\":0,\"areaSegments\":200,\"anchorLocalOffset\":{\"x\":0,\"y\":0},\"size\":{\"x\":1,\"y\":1}}],\"buildings\":["
#define TAIL "]}"

#define ASSEMBLER "{\"index\":%d,\"areaIndex\":0,\"localOffset\":[{\"x\":%lf,\"y\":%lf,\"z\":0},{\"x\":%lf,\"y\":%lf,\"z\":0}],\"yaw\":[0,0],\"itemId\":2305,\"modelIndex\":67,\"outputObjIdx\":-1,\"inputObjIdx\":-1,\"outputToSlot\":0,\"inputFromSlot\":0,\"outputFromSlot\":0,\"inputToSlot\":0,\"outputOffset\":0,\"inputOffset\":0,\"recipeId\":0,\"filterId\":0,\"parameters\":{\"acceleratorMode\":0}}"
#define SMELTER "{\"index\":%d,\"areaIndex\":0,\"localOffset\":[{\"x\":%lf,\"y\":%lf,\"z\":0},{\"x\":%lf,\"y\":%lf,\"z\":0}],\"yaw\":[0,0],\"itemId\":2315,\"modelIndex\":194,\"outputObjIdx\":-1,\"inputObjIdx\":-1,\"outputToSlot\":0,\"inputFromSlot\":0,\"outputFromSlot\":0,\"inputToSlot\":0,\"outputOffset\":0,\"inputOffset\":0,\"recipeId\":0,\"filterId\":0,\"parameters\":{\"acceleratorMode\":0}}"
#define CHEMICAL "{\"index\":%d,\"areaIndex\":0,\"localOffset\":[{\"x\":%lf,\"y\":%lf,\"z\":0},{\"x\":%lf,\"y\":%lf,\"z\":0}],\"yaw\":[0,0],\"itemId\":2317,\"modelIndex\":376,\"outputObjIdx\":-1,\"inputObjIdx\":-1,\"outputToSlot\":0,\"inputFromSlot\":0,\"outputFromSlot\":0,\"inputToSlot\":0,\"outputOffset\":0,\"inputOffset\":0,\"recipeId\":0,\"filterId\":0,\"parameters\":{\"acceleratorMode\":0}}"

// TODO 可自定义

#define Y_MAX 180.0
#define Y_MIN 0.0

// TODO 可自定义
#define SMELTER_X (2.31118)
#define SMELTER_Y (2.31118)
#define SMELTER2_X (2.31118)
#define SMELTER2_Y (5.31118)

#define ASSEMBLER_X (3.0416)
#define ASSEMBLER_Y (3.0416)
#define ASSEMBLER2_X (3.0416)
#define ASSEMBLER2_Y (6.7416)

#define CHEMICAL_X (6.6)
#define CHEMICAL_Y (3.6)
#define CHEMICAL2_X (6.6)
#define CHEMICAL2_Y (7.2)

#define CHROMOSOMES_LENGTH 64

typedef enum {
    none = 0,
    smelter2,
    assembler2,
    chemical2,
    ENUM_END
}module_t;

const double building_size_y[ENUM_END + 1] = {
    0.0,
    SMELTER_Y,
    ASSEMBLER_Y,
    CHEMICAL_Y,
    NAN
};

const double building_size_x[ENUM_END + 1] = {
    0.0,
    SMELTER_X,
    ASSEMBLER_X,
    CHEMICAL_X,
    NAN
};

const double module_size_y[ENUM_END + 1] = {
    0.0,
    SMELTER2_Y,
    ASSEMBLER2_Y,
    CHEMICAL2_Y,
    NAN
};

const double module_size_x[ENUM_END + 1] = {
    0.0,
    SMELTER2_X,
    ASSEMBLER2_X,
    CHEMICAL2_X,
    NAN
};

// 遗传算法的数据结构
typedef struct {
    module_t* gene;
    pcg32_random_t pcg;
    double score;
}individual_t;

// 获取建筑的对角线长度
double get_diagonal(double size_x, double size_y) {
    return sqrt(size_x * size_x + size_y * size_y);
}

// 获取蓝图坐标y处的半径
double get_r(double bp_y) {
    double radian_y = bp_y / 500.0 * M_PI;
    double r = cos(radian_y) * 500.0 / M_PI;
    return r;
}

// 输入建筑高纬度的中轴处对应的蓝图坐标y和建筑尺寸的xy，返回此处的蓝图最大x跨度
double get_blueprint_x_max(double bp_y, double size_x, double size_y) {
    // 先检查建筑对角线是否比当前纬度圆的半径还大，不然没有意义
    double module_diagonal = get_diagonal(size_x, size_y);
    if(fabs(bp_y) >= 250.0 - module_diagonal)
        return 0.0;
    double r = get_r(bp_y);
    return (2 * M_PI * r) / 10.0; // 1/20球
}

// 输入建筑低纬度的顶角处对应的蓝图坐标y和建筑尺寸的xy，返回此处建筑最大的y跨度
double get_module_y_max(double bp_y, double size_x, double size_y) {
    // 先检查建筑对角线是否比当前纬度圆的半径还大，不然没有意义
    double module_diagonal = get_diagonal(size_x, size_y);
    if(fabs(bp_y) >= 250.0 - module_diagonal)
        return module_diagonal;
    double r = get_r(bp_y);
    return r - sqrt(r * r - 0.25 * size_x * size_x) + size_y;
}

// 把基因转换成蓝图坐标
void get_bp_y(double bp_y[CHROMOSOMES_LENGTH], module_t gene[CHROMOSOMES_LENGTH]) {
    bp_y[0] = Y_MIN + get_module_y_max(0.0, module_size_x[gene[0]], module_size_y[gene[0]]);
    for(size_t i = 1; i < CHROMOSOMES_LENGTH; i++)
        bp_y[i] = bp_y[i - 1] + get_module_y_max(bp_y[i - 1], module_size_x[gene[i]], module_size_y[gene[i]]);
}

// 通过蓝图坐标和模型计算这行最多能放几个建筑
size_t get_module_count(double y, module_t id) {
    double bp_x_max = get_blueprint_x_max(y, module_size_x[id], module_size_y[id]);
    return (size_t)(bp_x_max / module_size_x[id]);
}

// 评价函数
double get_score(module_t gene[CHROMOSOMES_LENGTH], double* need) {
    // 先检查y是否超标，如果超标返回-y
    double bp_y[CHROMOSOMES_LENGTH] = { 0.0 };
    get_bp_y(bp_y, gene);
    if(bp_y[CHROMOSOMES_LENGTH - 1] > Y_MAX)
        return -bp_y[CHROMOSOMES_LENGTH - 1];

    // 否则计算建筑数，从赤道开始，向极地排列建筑
    size_t module_count[ENUM_END] = { 0 };
    for(size_t i = 0; i < CHROMOSOMES_LENGTH; i++) {
        if(gene[i] != none)
            module_count[gene[i]] += get_module_count(bp_y[i], gene[i]);
    }
    // 找出瓶颈的建筑是哪个
    double weight[ENUM_END] = { 0.0 };
    module_t bottleneck = 1;
    for(size_t i = 1; i < ENUM_END; i++)
        weight[i] = (double)module_count[i] / (double)need[i];
    for(size_t i = 2; i < ENUM_END; i++) {
        if(weight[i] < weight[bottleneck])
            bottleneck = i;
    }

    // 计算分数
    return weight[bottleneck];
}

void output_building(FILE* fp, const char* building, size_t* ptr_index, double x, double y) {
    if(*ptr_index > 0)
        fprintf(fp, ",");
    fprintf(fp, building, (*ptr_index)++, x, y, x, y);
}

void output_smelter2(FILE* fp, size_t* ptr_index, double x, double y) {
    double y_0 = y + module_size_y[smelter2] * 0.5 - building_size_y[smelter2] * 0.5;
    double y_1 = y - module_size_y[smelter2] * 0.5 + building_size_y[smelter2] * 0.5;
    output_building(fp, SMELTER, ptr_index, x, y_0);
    output_building(fp, SMELTER, ptr_index, x, y_1);
}

void output_assembler2(FILE* fp, size_t* ptr_index, double x, double y) {
    double y_0 = y + module_size_y[assembler2] * 0.5 - building_size_y[assembler2] * 0.5;
    double y_1 = y - module_size_y[assembler2] * 0.5 + building_size_y[assembler2] * 0.5;
    output_building(fp, ASSEMBLER, ptr_index, x, y_0);
    output_building(fp, ASSEMBLER, ptr_index, x, y_1);
}

void output_chemical2(FILE* fp, size_t* ptr_index, double x, double y) {
    double y_0 = y + module_size_y[chemical2] * 0.5 - building_size_y[chemical2] * 0.5 - 0.5;
    double y_1 = y - module_size_y[chemical2] * 0.5 + building_size_y[chemical2] * 0.5 - 0.5;
    output_building(fp, CHEMICAL, ptr_index, x, y_0);
    output_building(fp, CHEMICAL, ptr_index, x, y_1);
}

void output(module_t gene[CHROMOSOMES_LENGTH], double need[ENUM_END]) {
    // 输出基因的编码
    fprintf(stderr, "\ngenes:\n");
    for(size_t i = 0; i < CHROMOSOMES_LENGTH; i++) {
        if(gene[i] != none)
            fprintf(stderr, "%d", gene[i]);
    }
    double bp_y[CHROMOSOMES_LENGTH] = { 0.0 };
    get_bp_y(bp_y, gene);
    size_t module_count[ENUM_END] = { 0 };
    size_t module_count_gene[CHROMOSOMES_LENGTH] = { 0 };
    for(size_t i = 0; i < CHROMOSOMES_LENGTH; i++) {
        if(gene[i] != none) {
            module_count_gene[i] = get_module_count(bp_y[i], gene[i]);
            module_count[gene[i]] += module_count_gene[i];
        }
    }
    // 输出模块统计
    fprintf(stderr, "\nmodule count:\n");
    for(size_t i = 1; i < ENUM_END; i++)
        fprintf(stderr, "total %lld: %lld\t, score = %lf\n", i, module_count[i], (double)module_count[i] / need[i]);

    // 输出蓝图
    FILE* fp = stdout;
    fprintf(fp, HEAD);

    size_t index = 0;

    for(size_t i = 0; i < CHROMOSOMES_LENGTH; i++) {
        if(gene[i] != none) {
            double y_last = i > 0 ? bp_y[i - 1] : 0.0;
            double y_base = (bp_y[i] + y_last) * 0.5;
            for(size_t j = 0; j < module_count_gene[i]; j++) {
                double x_n = ((double)j + 0.5) / (double)module_count_gene[i] * 100.0; // 1/20球

                switch(gene[i]) {
                case smelter2:
                output_smelter2(fp, &index, x_n, y_base);
                break;

                case  assembler2:
                output_assembler2(fp, &index, x_n, y_base);
                break;

                case chemical2:
                output_chemical2(fp, &index, x_n, y_base);
                break;

                default:
                fprintf(stderr, "!!! Warning: unknow module !!!\n");
                }
            }
        }
    }

    fprintf(fp, TAIL);
}

// 遗传算法的种群生成
void init_individual(individual_t* ptr, double* need, uint64_t seed) {
    // 初始化伪随机数生成器
    pcg32_random_t pcg = { seed, 0 };
    ptr->pcg = pcg;
    pcg32_random_r(&(ptr->pcg));
    // 对need求和
    need[0] = 0.0;
    for(size_t i = 1; i < ENUM_END; i++)
        need[0] += need[i];
    // 初始化初代个体
    ptr->gene = calloc(CHROMOSOMES_LENGTH, sizeof(module_t));
    for(size_t i = 0; i < CHROMOSOMES_LENGTH; i++)
        ptr->gene[i] = pcg32_random_r(&(ptr->pcg)) % ENUM_END;
    // 计算个体分数
    ptr->score = get_score(ptr->gene, need);
}

// 遗传算法的进化过程

int cmp(const void* ptr_a, const void* ptr_b) {
    const individual_t* a = (const individual_t*)ptr_a;
    const individual_t* b = (const individual_t*)ptr_b;
    if(a->score < b->score)
        return -1;
    if(a->score > b->score)
        return 1;
    return 0;
}

int main(int argc, char* argv[]) {
    const time_t time_start = time(NULL);

    double need[ENUM_END] = { 0, 43, 300, 57 }; // 第一个数必须填0

    // 种群生成
    size_t individual_count = 65536;
    individual_t* individual = calloc(individual_count, sizeof(individual_t));
    for(size_t i = 0; i < individual_count; i++)
        init_individual(&(individual[i]), need, (uint64_t)((size_t)time_start + i + 1));

    // 适应度排序，适应度越高越靠后
    qsort(individual, individual_count, sizeof(individual_t), cmp);

    // 迭代进化
    const double p_crossover = 0.5;
    const double p_mutation = 0.02;
    pcg32_random_t pcg = { (uint64_t)time_start, 0 };
    pcg32_random_r(&pcg);

    size_t kill_count = (size_t)((double)individual_count * p_crossover);
    for(size_t iteration = 0; iteration < 256; iteration++) {
        // 交叉
    #pragma omp parallel for
        for(size_t i = 0; i < kill_count; i++) {
            size_t mom = kill_count + pcg32_random_r(&pcg) % (individual_count - kill_count);
            size_t dad = kill_count + pcg32_random_r(&pcg) % (individual_count - kill_count);
            size_t crossover_point = pcg32_random_r(&pcg) % (CHROMOSOMES_LENGTH - 1) + 1;
            for(size_t j = 0; j < crossover_point; j++)
                individual[i].gene[j] = individual[mom].gene[j];
            for(size_t j = crossover_point; j < CHROMOSOMES_LENGTH; j++)
                individual[i].gene[j] = individual[dad].gene[j];
            individual[i].score = get_score(individual[i].gene, need);
        }
        // 变异
    #pragma omp parallel for
        for(size_t i = 0; i < individual_count; i++) {
            if(((double)pcg32_random_r(&pcg) / (double)UINT32_MAX) < p_mutation) {
                size_t mutation_point = pcg32_random_r(&pcg) % CHROMOSOMES_LENGTH;
                individual[i].gene[mutation_point] = pcg32_random_r(&pcg) % ENUM_END;
                individual[i].score = get_score(individual[i].gene, need);
            }
        }
        // 适应度排序，适应度越高越靠后
        qsort(individual, individual_count, sizeof(individual_t), cmp);
        fprintf(stderr, "iteration = %lld, best = %lf\n", iteration, individual[individual_count - 1].score);
    }

    // 输出
    output(individual[individual_count - 1].gene, need);

    // free
    for(size_t i = 0; i < individual_count; i++)
        free(individual->gene);
    free(individual);
    return 0;
}