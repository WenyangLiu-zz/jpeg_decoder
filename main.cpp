#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <map>
#include <math.h>
#include "qdbmp.h"
#include <vector>
#include <iostream>
#include <time.h>
#include <string>
#include <bitset>
#include <iomanip> 

using namespace std;

const int SOI_MARKER = 0xD8;
const int APP0_MARKER = 0xE0;
const int DQT_MARKER = 0xDB;
const int SOF_MARKER = 0xC0;
const int DHT_MARKER = 0xC4;
const int SOS_MARKER = 0xDA;
const int EOI_MARKER = 0xD9;
const int COM_MARKER = 0xFE;
const int DRI_MARKER = 0xDD;

struct {
    int height;
    int width;
} image;

struct {
    unsigned char id;
    unsigned char width;
    unsigned char height;
    unsigned char quant;
} subVector[4];

unsigned char maxWidth, maxHeight;

struct acCode {
    unsigned char len;    /*code length*/
    unsigned char zeros;  /*zero run*/
    int value;
};

struct RGB {
    unsigned char R, G, B;
};

typedef double BLOCK[8][8];

int quantTable[4][128];

const int DC = 0;
const int AC = 1;
map<pair<unsigned char, unsigned int>, unsigned char> huffTable[2][2];
map<int, pair<int, string>> byteStream;
double cos_cache[200];
int stream_num = 1;
int dc[4] = {0, 0, 0, 0};  


void init_cos_cache() {
    for (int i = 0; i < 200; i++) {
        cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

class MCU {
    public:
    BLOCK mcu[4][2][2];
    // 除錯用
    void show() {
        printf("*************** mcu show ***********************\n");
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    printf("mcu id: %d, %d %d\n", id, h, w);
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            printf("%lf ", mcu[id][h][w][i][j]);
                        }
                    }
                }
            }
        }
    };
    void quantify() {
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] *= quantTable[subVector[id].quant][i*8 + j];
                        }
                    }
                }
            }
        }
    };
    void zigzag() {
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    int zz[8][8] = {
                            { 0,  1,  5,  6, 14, 15, 27, 28},
                            { 2,  4,  7, 13, 16, 26, 29, 42},
                            { 3,  8, 12, 17, 25, 30, 41, 43},
                            { 9, 11, 18, 24, 31, 40, 44, 53},
                            {10, 19, 23, 32, 39, 45, 52, 54},
                            {20, 22, 33, 38, 46, 51, 55, 60},
                            {21, 34, 37, 47, 50, 56, 59, 61},
                            {35, 36, 48, 49, 57, 58, 62, 63}
                    };
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            zz[i][j] = mcu[id][h][w][zz[i][j] / 8][zz[i][j] % 8];
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] = zz[i][j];
                        }
                    }
                }
            }
        }
    };
    void idct() {
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    double tmp[8][8] = {0};
//					 照定義展開，效能低下
//                     for (int i = 0; i < 8; i++) {
//                         for (int j = 0; j < 8; j++) {
//                             for (int x = 0; x < 8; x++) {
//                                 for (int y = 0; y < 8; y++) {
//                                     tmp[i][j] += (cc(x, y) * mcu[id][h][w][x][y] * cos((2*i+1)*M_PI/16.0*x) * cos((2*j+1)*M_PI/16.0*y));
//                                 }
//                             }
//                             tmp[i][j] /= 4.0;
//                         }
//                     }
                    // 計算兩次一維idct去計算二維idct
                    double s[8][8] = {};
                    for (int j = 0; j < 8; j++) {
                        for (int x = 0; x < 8; x++) {
                            for (int y = 0; y < 8; y++) {
                                s[j][x] += c (y) * mcu[id][h][w][x][y] * cos_cache[(j + j + 1) * y];
                            }
                            s[j][x] = s[j][x] / 2.0;
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            for (int x = 0; x < 8; x++) {
                                tmp[i][j] += c(x) * s[j][x] * cos_cache[(i + i + 1) * x];
                            }
                            tmp[i][j] = tmp[i][j] / 2.0;
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] = tmp[i][j];
                        }
                    }
                }
            }
        }
    }
    void decode() {
        this->quantify();
        this->zigzag();
        this->idct();
    }
    RGB **toRGB() {
        RGB **ret = (RGB **)malloc(sizeof(RGB **) * maxHeight * 8);
        for (int i = 0; i < maxHeight * 8; i++) {
            ret[i] = (RGB *)malloc(sizeof(RGB *) * maxWidth * 8);
        }
        for (int i = 0; i < maxHeight * 8; i++) {
            for (int j = 0; j < maxWidth * 8; j++) {
                double Y = trans(1, i, j);
                double Cb = trans(2, i, j);
                double Cr = trans(3, i, j);
                ret[i][j].R = chomp(Y + 1.402*Cr + 128);
                ret[i][j].G = chomp(Y - 0.34414*Cb - 0.71414*Cr + 128);
                ret[i][j].B = chomp(Y + 1.772*Cb + 128);
            }
        }
        return ret;
    }
private:
    double cc(int i, int j) {
        if (i == 0 && j == 0) {
            return 1.0/2.0;
        } else if (i == 0 || j == 0) {
            return 1.0/sqrt(2.0);
        } else {
            return 1.0;
        }
    }
    double c(int i) {
        static double x = 1.0/sqrt(2.0);
        if (i == 0) {
            return x;
        } else {
            return 1.0;
        }
    }
    unsigned char chomp(double x) {
        if (x > 255.0) {
            return 255;
        } else if (x < 0) {
            return 0;
        } else {
            return (unsigned char) x;
        }
    }
    double trans(int id, int h, int w) {
        int vh = h * subVector[id].height / maxHeight;
        int vw = w * subVector[id].width / maxWidth;
        return mcu[id][vh / 8][vw / 8][vh % 8][vw % 8];
    }
};

// 讀取 Section 用之輔助函式

void showSectionName(const char *s) {
    printf("************************ %s **************************\n", s);
    return;
}



void outputMap() { //迭代器里数据除去了一些多余的0xFF
    map<int, pair<int, string>>::iterator iter;
    int i = 0;
	for (iter = byteStream.begin(); iter != byteStream.end(); iter++) {
        cout << "Start:" <<iter->first << endl;
        cout << "Cur bit position:" << iter->second.first << endl;
        for(unsigned char x : iter->second.second) {
            if(i % 32 == 0) {
                cout << endl;
            }
            cout << " 0x" << hex << setfill('0') << setw(2) << (int)x;
            i++;
        }  
	}
    iter = byteStream.begin();
    int start = iter->first;
    cout << endl << dec << iter->second.second.size() << " " << start ;
    for (int i = 0;i< 2; i++) {
        unsigned char x = iter->second.second[50000+i];
        cout << " 0x" << hex << setfill('0') << setw(2) << (int)x;
    }
    cout << endl;
}

unsigned int readSectionLength(FILE *f) {
    unsigned char c;
    unsigned int length;
    fread(&c, 1, 1, f);
    length = c;
    fread(&c, 1, 1, f);
    length = length * 256 + c;
    return length;
}

unsigned int EnterNewSection(FILE *f, const char *s) {
    showSectionName(s);
    unsigned int len = readSectionLength(f);
    printf("本區段長度為 %d\n", len);
    return len;
}

// 讀取各 Section 之函式

void readCOM(FILE *f) {
    unsigned int len = EnterNewSection(f, "COM");
    unsigned char c;
    for (int i = 0; i < len - 2; i++) {
        fread(&c, 1, 1, f);
        printf("%c", c);
    }
    printf("\n");
}

void readAPP(FILE *f) {
    unsigned int len = EnterNewSection(f, "APP0");
    char m[5];
    fread(m, 1, 5, f);
    printf("使用 %s\n", m);
    unsigned char v[2];
    fread(v, 1, 2, f);
    printf("版本 %d.%d\n", v[0], v[1]);
    fseek(f, 1, SEEK_CUR);
    fread(v, 1, 2, f);
    printf("x方向像素密度：%d\n", v[0] * 16 + v[1]);
    fread(v, 1, 2, f);
    printf("y方向像素密度：%d\n", v[0] * 16 + v[1]);
    fseek(f, len - 14, SEEK_CUR);
}
void readDQT(FILE *f) {
    unsigned int len = EnterNewSection(f, "DQT");
    len -= 2;
    while (len > 0) {
        unsigned char c;
        fread(&c, 1, 1, f);
        len--;
        unsigned precision = c >> 4 == 0 ? 8 : 16;
        printf("精度：%d\n", precision);
        precision /= 8;
        unsigned char id = c & 0x0F;
        printf("量化表ID: %d\n", id);
        for (int i = 0; i < 64; i++) {
            unsigned char t = 0;
            for (int p = 0; p < precision; p++) {
                unsigned char s;
                fread(&s, 1, 1, f);
                t == t << 8;
                t += s;
            }
            quantTable[id][i] = t;
        }
        for (int i = 0; i < 64; i++) {
            if (i % 8 == 0) {
                printf("\n");
            }
            printf("%2d ", quantTable[id][i]);
        }
        printf("\n");
        len -= (precision*64);
    }
}

void readSOF(FILE *f) {
    unsigned int len = EnterNewSection(f, "SOF");
    fseek(f, 1, SEEK_CUR); // 精度
    unsigned char v[3];
    fread(v, 1, 2, f);
    // TODO: 高度跟寬度不確定
    image.height = v[0] * 256 + v[1];
    fread(v, 1, 2, f);
    image.width = v[0] * 256 + v[1];
    printf("高*寬: %d*%d\n", image.height, image.width);
    fseek(f, 1, SEEK_CUR); // 顏色分量數，固定為3
    for (int i = 0; i < 3; i++) {
        fread(v, 1, 3, f);
        printf("顏色分量ID：%d\n", v[0]);
        printf("水平採樣因子：%d\n", v[1] >> 4);
        printf("垂直採樣因子：%d\n", v[1] & 0x0F);
        printf("量化表ID：%d\n", v[2]);
        subVector[v[0]].id = v[0];
        subVector[v[0]].width = v[1] >> 4;
        subVector[v[0]].height = v[1] & 0x0F;
        subVector[v[0]].quant = v[2];
        maxHeight = (maxHeight > subVector[v[0]].height ? maxHeight : subVector[v[0]].height);
        maxWidth = (maxWidth > subVector[v[0]].width ? maxWidth : subVector[v[0]].width);
    }
}

//创建范式huffcode，与标准的huffmancode不一样的点在与，加了3点约束
//1. 最小编码长度的第一个编码必须从 0 开始。
//2. 相同长度编码必须是连续的。
//3. 编码长度为j的第一个符号可以从编码长度为j-1的最后一个符号所得知，即c_j = 2(c_{j-1}+1) -> c_j = (c_{j-1}+1) << 1 加一左移

pair<unsigned char, unsigned int>* createHuffCode(unsigned char *a, unsigned int number) {
    int si = sizeof(pair<unsigned char, unsigned int>);
    auto ret = (pair<unsigned char, unsigned int>*)malloc(si*number);        //分配一定数目的空间用来存放huffcode到ret中
    int code = 0;
    int count = 0;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < a[i]; j++) {                //数组a中存放了从1bit到16bit中，码字的个数
            ret[count++] = make_pair(i + 1, code);
            code += 1;               //第二条约束
        }
        code = code << 1;                //第三条约束      
    }
    return ret;
}

void readDHT(FILE *f) {
    unsigned int len = EnterNewSection(f, "DHT");
    len -= 2;
    while (len > 0) {
        unsigned char v;      //接下来的读取1个字节，高4位表示AC或者DC哈夫曼表。低4位表示哈夫曼表的ID，这里ID要么0要么1。
        fread(&v, 1, 1, f);              //int 4个字节，char 是一个字节， char i = 1 -> 0000 0001
        unsigned char DCorAC = v >> 4;       
        printf(DCorAC == 0 ? "DC " : "AC ");
        unsigned char id = v & 0x0F;
        printf("%d\n", id);

        unsigned char a[16];
        fread(a, 1, 16, f);
        unsigned int number = 0;
        for (int i = 0; i < 16; i++) {
            printf("%d ", a[i]);
            number += a[i];
        }
        printf("不同位宽下的数目\n");
        auto huffCode = createHuffCode(a, number);             //huffcode 码表，其中形式就是pair数组，每个pair是一个<char, int>
        for (int i = 0; i < number; i++) {                     //number指的是码表里面的个数，每次取一个元素放到码表里
            unsigned char v;
            fread(&v, 1, 1, f);                                //读取1个byte
            huffTable[DCorAC][id][huffCode[i]] = v;
            //printf("%d %d: %d\n", huffCode[i].first, huffCode[i].second, v);           //这一步是将码表打印出来 -------- Code length  | Number | Symbol  Symbol是保存下一次取出Code 的位数
            bitset<16> codeword(huffCode[i].second);
            cout << "(" << std::dec << (int)huffCode[i].first << "," << codeword.to_string() << "): 0x" << hex << (int)v << endl;
        }
        free(huffCode); 
        len -= (1 + 16 + number);
    }
}
void readSOS(FILE *f) {
    unsigned int len = EnterNewSection(f, "SOS");
    fseek(f, 1, SEEK_CUR);   // 顏色分量數，固定為3
    for (int i = 0; i < 3; i++) {
        unsigned char v[1];
        fread(v, 1, 1, f);
        printf("顏色分量id：%d\n", v[0]);
        fread(v, 1, 1, f);
        printf("DC霍夫曼id：%d\n", v[0] >> 4);
        printf("AC霍夫曼id：%d\n", v[0] & 0x0F);
    }
    fseek(f, 3, SEEK_CUR);
}


// 必須連續呼叫getBit，中間被fread斷掉就會出問題，每次读取一个bit， byteStream[key].first 表示的当前的bit值
bool getBit(int key) {
    // 虽然每次读取一个字节，但是利用了count的循环，控制了一个字节可以获取8个bit
    int cur_byte = byteStream[key].first / 8;
    int cur_bit = byteStream[key].first % 8;
    string bit_string = byteStream[key].second;
    byteStream[key].first += 1;
    return byteStream[key].second[cur_byte] & (1 << (7 - cur_bit));
}

void removeBit(int key, int num){
    for(int i=0;i<num; i++){
        getBit(key);
    }
}

unsigned char matchHuff(unsigned char number, unsigned char ACorDC, int key) {
    char codeword[20] = "";
    unsigned int codeword_val = 0;
    unsigned char codeLen;
    for (int count = 1; ; count++) {  
        bool cur_bit = getBit(key);
        codeword_val = codeword_val << 1;
        codeword_val += (unsigned int)cur_bit;     //每次读取一个bit，每次左移，可以得到其中
        codeword[count - 1] = '0' + (unsigned int)cur_bit;
        // 迭代器找到了该元素，没找到会返回end迭代器
        if (huffTable[ACorDC][number].find(std::make_pair(count, codeword_val)) != huffTable[ACorDC][number].end()) { 
            codeLen = huffTable[ACorDC][number][std::make_pair(count, codeword_val)];
            return codeLen;
        }

        if (ACorDC == DC) {
            // DC codeword 最大为11个元素，如果失败了，缺失key，重来。此时文件头已经读过去了,抛弃11bit
            if (count > 6) { //该图像是6bit
                codeword[11] = '\0';
                printf("Byte position %d %d %d %s %d\n", key, byteStream[key].first / 8, byteStream[key].first % 8, codeword, ACorDC);
                fprintf(stderr, "DC key not found\n"); 
                count = 1; codeword_val = 0; 

            }
        } else {
            // AC codeword 最大为16个元素，如果失败了，缺失key，重来。此时文件头已经读过去了，抛弃了16个bits
            if (count > 16) {
                codeword[15] = '\0';
                printf("Byte position %d %d %d %s %d\n", key, byteStream[key].first / 8, byteStream[key].first % 8, codeword, ACorDC);
                fprintf(stderr, "AC key not found\n"); 
                count = 1; codeword_val = 0;
            }
        }
    }
}

int readDC(unsigned char number, int key) {
    unsigned char codeLen = matchHuff(number, DC, key);  //查表,得到codelen，即下一次取多少bit
    if (codeLen == 0) { return 0; }  
    unsigned char first = getBit(key); //符号位
    int ret = 1;
    for (int i = 1; i < codeLen; i++) {
        unsigned char b = getBit(key);
        ret = ret << 1;
        ret += first ? b : !b;
    }
    ret = first ? ret : -ret;
    return ret;
}

// 計算ZRL
acCode readAC(unsigned char number, int key) {
    unsigned char x = matchHuff(number, AC, key); //查表，返回1 byte
    unsigned char zeros = x >> 4;               //前面的0
    unsigned char codeLen = x & 0x0F;           //编码的长度
    if (x == 0) {   //EOB
        return acCode{0,0,0};  
    } else if (x == 0xF0) {     // 连续16个0
        return acCode{0, 16, 0};
    }
    unsigned char first = getBit(key);
    int value = 1;
    /* value = 2^(codelen - 1) + signed offset = 4  右移（codelen-1） + 剩余部分
       1100
       1 - true
       1 - 右移动加1 - 11
       0 - 110
       0 - 1100 = 14
    */
    for (int i = 1; i < codeLen; i++) {
        unsigned char b = getBit(key);
        value = value << 1;       // 右移
        value += first ? b : !b;  // 如果first为真，就取b，否则！b         
    }
    value = first ? value : -value;
    return acCode{codeLen, zeros, value};
}

MCU readMCU(int key) {
    auto mcu = MCU();
    for (int i = 1; i <= 3; i++) {
        for (int h = 0; h < subVector[i].height; h++) {
            for (int w = 0; w < subVector[i].width; w++) {
                int DCvalue = readDC(i/2, key);
                dc[i] = DCvalue + dc[i];
                // if(i == 1) {
                //     FILE *fp = fopen("o.txt", "a");
                //     fprintf(fp,"%d\n",DCvalue);  //字符使用%c
                //     fclose(fp);
                // }
                mcu.mcu[i][h][w][0][0] = dc[i];
                unsigned int count = 1;
                while (count < 64) {
                    acCode ac = readAC(i/2, key);
                    /*zeros最多只可以有16个0*/
                    if (ac.len == 0 && ac.zeros == 16) {
                        for (int j = 0; j < ac.zeros; j++) {
                            mcu.mcu[i][h][w][count/8][count%8] = 0;
                            count++;
                        }
                    } else if (ac.len == 0 && ac.zeros == 0) {
                        /*如果ac碰到了EOB，则给剩下的值赋予0*/
                        while (count < 64) {
                            mcu.mcu[i][h][w][count/8][count%8] = 0;
                            count++;
                        }
                        break;
                    } else {
                        // 对前面的0元素赋0，以及当前元素赋予初值value
                        for (int j = 0; j < ac.zeros; j++) {
                            mcu.mcu[i][h][w][count/8][count%8] = 0;
                            count++;
                        }
                        mcu.mcu[i][h][w][count/8][count%8] = ac.value;  //给mcu块赋值
                        count++;
                    }
                }
            
            }
        }
    }
    return mcu;
}

void readData() {
    printf("Read data \n");
    int w = (image.width - 1) / (8*maxWidth) + 1;
    int h = (image.height - 1) / (8*maxHeight) + 1;
    int inter_mcu = w * h / stream_num; 
    int index = 0;
    BMP *bmp = BMP_Create(maxWidth * 8 * w, maxHeight * 8 * h, 24);
    map<int, pair<int, string>>::iterator iter = byteStream.begin();
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            MCU mcu = readMCU(iter->first);
            mcu.decode();
            RGB **b = mcu.toRGB();
            for (int y = i*8*maxHeight; y < (i+1)*8*maxHeight; y++) {
                for (int x = j*8*maxWidth; x < (j+1)*8*maxWidth; x++) {
                    int by = y - i*8*maxHeight;
                    int bx = x - j*8*maxWidth;
                    BMP_SetPixelRGB(bmp, x, y, b[by][bx].R, b[by][bx].G, b[by][bx].B);
                }
            }
            int mcu_count = i*w + j+ 1;
            if(mcu_count % inter_mcu == 0) {
                iter ++; index ++;
                dc[0] = 0; dc[1] = 0; dc[2]= 0;  dc[3]= 0; 
                cout << index << " " << stream_num << endl;
            }
        }
    }
    BMP_WriteFile(bmp, "out.bmp");
    outputMap();
}


unsigned int readDRI(FILE *f) {
    unsigned int len = EnterNewSection(f, "DRI");
    unsigned char c[2];
    fread(&c, 1, 2, f);
    unsigned int interv  = c[0] << 8 + c[1];
    printf("复位的间隔是 %d\n", interv);
    return interv;
}


void matchRst(FILE *f){
    printf("scan for the image data\n");
    bool flag = true; 
    string stream;
    int start_pos_key = ftell(f);
    unsigned char cur;
    while (flag) {
        fread(&cur, 1, 1, f);
        if (cur != 0xFF) {
            stream.push_back(cur);
            continue;
        }
        // 到这里时 curr == 0xFF
        // 下一个字节为 next == 0x00, 0xFF, 0xab
        unsigned char next;
        fread(&next, 1, 1, f);
        while(next== 0xFF) { // 当check不为FF跳出
            fread(&next, 1, 1, f);
        }   

        if(next == 0x00) {
            stream.push_back(0xFF);
            continue;
        }
        // 把想要的元素放到map里
        std::pair<int, string> val = std::make_pair(0, stream);
        byteStream[start_pos_key] = val;
        stream.clear(); // 清除所有的元素
        // 能到这里说明是碰到了标记位 0xFFab 
        switch (next)  
        {
        case 0xD0:
        case 0xD1:
        case 0xD2:
        case 0xD3:
        case 0xD4:
        case 0xD5:
        case 0xD6:
        case 0xD7: // 如果碰到rst标记
            flag = true;
            start_pos_key = ftell(f); //更新key
            break;
        default: //其他标记，直接跳出
            if(next != 0xD9){
                fprintf(stderr, "data段有不是0xFF00的数据 %02x\n", next);
            }
            fseek(f, -2, SEEK_CUR);
            flag = false;
        }
    } 
    //readData(f);
    stream_num = byteStream.size();
}

void readStream(FILE *f) {
    int lock = 1;
    int interv = 0;
    unsigned char c;
    while (lock) {
        fread(&c, 1, 1, f);
        if(c != 0xFF) {
            continue;
        }
        fread(&c, 1, 1, f);
        switch (c) {
            case SOI_MARKER:
                printf("Start of Image\n");
                break;
            case APP0_MARKER:
            case 0xE1:
            case 0xE2:
            case 0xE3:
            case 0xE4:
            case 0xE5:
            case 0xE6:
            case 0xE7:
            case 0xE8:
            case 0xE9:
            case 0xEA:
            case 0xEB:
            case 0xEC:
            case 0xED:
            case 0xEE:
            case 0xEF:
                readAPP(f);
                break;
            case COM_MARKER:
                readCOM(f);
                break;
            case DQT_MARKER:
                readDQT(f);
                break;
            case SOF_MARKER:
                readSOF(f);
                break;
            case DHT_MARKER:
                readDHT(f);
                break;
            case SOS_MARKER:
                readSOS(f);
                matchRst(f);
                //outputMap();
                readData();
                break;
            case DRI_MARKER:
                interv = readDRI(f);
                printf("%d --- \n", interv);
                break;
            case EOI_MARKER:
                lock = 0;
                printf("End of Image\n");
                break;
            default:
                printf("other marker: %02x\n", c);
        }
    }
    if (fread(&c, 1, 1, f) != 0) {
        fprintf(stderr, "沒有吃完就結束\n");
    }
}

int main(int argc, char *argv[]) {
    clock_t start,end;     //定义clock_t变量
    start = clock();       //开始时间
    if (argc != 2) {
        fprintf(stderr, "用法：jpeg_decoder <jpeg file>\n");
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "檔案開啟失敗\n");
    }
    init_cos_cache();
    readStream(f);
    end = clock();
    printf("Time is %fs\n", double(end-start)/CLOCKS_PER_SEC);
    return 0;
}
