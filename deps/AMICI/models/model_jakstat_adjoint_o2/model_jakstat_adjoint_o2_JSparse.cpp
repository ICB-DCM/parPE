
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint_o2{

void JSparse_model_jakstat_adjoint_o2(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 9;
  JSparse->indexvals[3] = 10;
  JSparse->indexvals[4] = 54;
  JSparse->indexvals[5] = 55;
  JSparse->indexvals[6] = 63;
  JSparse->indexvals[7] = 64;
  JSparse->indexvals[8] = 72;
  JSparse->indexvals[9] = 73;
  JSparse->indexvals[10] = 81;
  JSparse->indexvals[11] = 82;
  JSparse->indexvals[12] = 90;
  JSparse->indexvals[13] = 91;
  JSparse->indexvals[14] = 1;
  JSparse->indexvals[15] = 2;
  JSparse->indexvals[16] = 10;
  JSparse->indexvals[17] = 11;
  JSparse->indexvals[18] = 19;
  JSparse->indexvals[19] = 20;
  JSparse->indexvals[20] = 28;
  JSparse->indexvals[21] = 29;
  JSparse->indexvals[22] = 37;
  JSparse->indexvals[23] = 38;
  JSparse->indexvals[24] = 46;
  JSparse->indexvals[25] = 47;
  JSparse->indexvals[26] = 55;
  JSparse->indexvals[27] = 56;
  JSparse->indexvals[28] = 64;
  JSparse->indexvals[29] = 65;
  JSparse->indexvals[30] = 73;
  JSparse->indexvals[31] = 74;
  JSparse->indexvals[32] = 82;
  JSparse->indexvals[33] = 83;
  JSparse->indexvals[34] = 91;
  JSparse->indexvals[35] = 92;
  JSparse->indexvals[36] = 100;
  JSparse->indexvals[37] = 101;
  JSparse->indexvals[38] = 109;
  JSparse->indexvals[39] = 110;
  JSparse->indexvals[40] = 118;
  JSparse->indexvals[41] = 119;
  JSparse->indexvals[42] = 127;
  JSparse->indexvals[43] = 128;
  JSparse->indexvals[44] = 136;
  JSparse->indexvals[45] = 137;
  JSparse->indexvals[46] = 145;
  JSparse->indexvals[47] = 146;
  JSparse->indexvals[48] = 154;
  JSparse->indexvals[49] = 155;
  JSparse->indexvals[50] = 2;
  JSparse->indexvals[51] = 3;
  JSparse->indexvals[52] = 29;
  JSparse->indexvals[53] = 30;
  JSparse->indexvals[54] = 3;
  JSparse->indexvals[55] = 4;
  JSparse->indexvals[56] = 39;
  JSparse->indexvals[57] = 40;
  JSparse->indexvals[58] = 4;
  JSparse->indexvals[59] = 5;
  JSparse->indexvals[60] = 40;
  JSparse->indexvals[61] = 41;
  JSparse->indexvals[62] = 5;
  JSparse->indexvals[63] = 6;
  JSparse->indexvals[64] = 41;
  JSparse->indexvals[65] = 42;
  JSparse->indexvals[66] = 6;
  JSparse->indexvals[67] = 7;
  JSparse->indexvals[68] = 42;
  JSparse->indexvals[69] = 43;
  JSparse->indexvals[70] = 7;
  JSparse->indexvals[71] = 8;
  JSparse->indexvals[72] = 43;
  JSparse->indexvals[73] = 44;
  JSparse->indexvals[74] = 0;
  JSparse->indexvals[75] = 8;
  JSparse->indexvals[76] = 36;
  JSparse->indexvals[77] = 44;
  JSparse->indexvals[78] = 9;
  JSparse->indexvals[79] = 10;
  JSparse->indexvals[80] = 10;
  JSparse->indexvals[81] = 11;
  JSparse->indexvals[82] = 11;
  JSparse->indexvals[83] = 12;
  JSparse->indexvals[84] = 12;
  JSparse->indexvals[85] = 13;
  JSparse->indexvals[86] = 13;
  JSparse->indexvals[87] = 14;
  JSparse->indexvals[88] = 14;
  JSparse->indexvals[89] = 15;
  JSparse->indexvals[90] = 15;
  JSparse->indexvals[91] = 16;
  JSparse->indexvals[92] = 16;
  JSparse->indexvals[93] = 17;
  JSparse->indexvals[94] = 9;
  JSparse->indexvals[95] = 17;
  JSparse->indexvals[96] = 18;
  JSparse->indexvals[97] = 19;
  JSparse->indexvals[98] = 19;
  JSparse->indexvals[99] = 20;
  JSparse->indexvals[100] = 20;
  JSparse->indexvals[101] = 21;
  JSparse->indexvals[102] = 21;
  JSparse->indexvals[103] = 22;
  JSparse->indexvals[104] = 22;
  JSparse->indexvals[105] = 23;
  JSparse->indexvals[106] = 23;
  JSparse->indexvals[107] = 24;
  JSparse->indexvals[108] = 24;
  JSparse->indexvals[109] = 25;
  JSparse->indexvals[110] = 25;
  JSparse->indexvals[111] = 26;
  JSparse->indexvals[112] = 18;
  JSparse->indexvals[113] = 26;
  JSparse->indexvals[114] = 27;
  JSparse->indexvals[115] = 28;
  JSparse->indexvals[116] = 28;
  JSparse->indexvals[117] = 29;
  JSparse->indexvals[118] = 29;
  JSparse->indexvals[119] = 30;
  JSparse->indexvals[120] = 30;
  JSparse->indexvals[121] = 31;
  JSparse->indexvals[122] = 31;
  JSparse->indexvals[123] = 32;
  JSparse->indexvals[124] = 32;
  JSparse->indexvals[125] = 33;
  JSparse->indexvals[126] = 33;
  JSparse->indexvals[127] = 34;
  JSparse->indexvals[128] = 34;
  JSparse->indexvals[129] = 35;
  JSparse->indexvals[130] = 27;
  JSparse->indexvals[131] = 35;
  JSparse->indexvals[132] = 36;
  JSparse->indexvals[133] = 37;
  JSparse->indexvals[134] = 37;
  JSparse->indexvals[135] = 38;
  JSparse->indexvals[136] = 38;
  JSparse->indexvals[137] = 39;
  JSparse->indexvals[138] = 39;
  JSparse->indexvals[139] = 40;
  JSparse->indexvals[140] = 40;
  JSparse->indexvals[141] = 41;
  JSparse->indexvals[142] = 41;
  JSparse->indexvals[143] = 42;
  JSparse->indexvals[144] = 42;
  JSparse->indexvals[145] = 43;
  JSparse->indexvals[146] = 43;
  JSparse->indexvals[147] = 44;
  JSparse->indexvals[148] = 36;
  JSparse->indexvals[149] = 44;
  JSparse->indexvals[150] = 45;
  JSparse->indexvals[151] = 46;
  JSparse->indexvals[152] = 46;
  JSparse->indexvals[153] = 47;
  JSparse->indexvals[154] = 47;
  JSparse->indexvals[155] = 48;
  JSparse->indexvals[156] = 48;
  JSparse->indexvals[157] = 49;
  JSparse->indexvals[158] = 49;
  JSparse->indexvals[159] = 50;
  JSparse->indexvals[160] = 50;
  JSparse->indexvals[161] = 51;
  JSparse->indexvals[162] = 51;
  JSparse->indexvals[163] = 52;
  JSparse->indexvals[164] = 52;
  JSparse->indexvals[165] = 53;
  JSparse->indexvals[166] = 45;
  JSparse->indexvals[167] = 53;
  JSparse->indexvals[168] = 54;
  JSparse->indexvals[169] = 55;
  JSparse->indexvals[170] = 55;
  JSparse->indexvals[171] = 56;
  JSparse->indexvals[172] = 56;
  JSparse->indexvals[173] = 57;
  JSparse->indexvals[174] = 57;
  JSparse->indexvals[175] = 58;
  JSparse->indexvals[176] = 58;
  JSparse->indexvals[177] = 59;
  JSparse->indexvals[178] = 59;
  JSparse->indexvals[179] = 60;
  JSparse->indexvals[180] = 60;
  JSparse->indexvals[181] = 61;
  JSparse->indexvals[182] = 61;
  JSparse->indexvals[183] = 62;
  JSparse->indexvals[184] = 54;
  JSparse->indexvals[185] = 62;
  JSparse->indexvals[186] = 63;
  JSparse->indexvals[187] = 64;
  JSparse->indexvals[188] = 64;
  JSparse->indexvals[189] = 65;
  JSparse->indexvals[190] = 65;
  JSparse->indexvals[191] = 66;
  JSparse->indexvals[192] = 66;
  JSparse->indexvals[193] = 67;
  JSparse->indexvals[194] = 67;
  JSparse->indexvals[195] = 68;
  JSparse->indexvals[196] = 68;
  JSparse->indexvals[197] = 69;
  JSparse->indexvals[198] = 69;
  JSparse->indexvals[199] = 70;
  JSparse->indexvals[200] = 70;
  JSparse->indexvals[201] = 71;
  JSparse->indexvals[202] = 63;
  JSparse->indexvals[203] = 71;
  JSparse->indexvals[204] = 72;
  JSparse->indexvals[205] = 73;
  JSparse->indexvals[206] = 73;
  JSparse->indexvals[207] = 74;
  JSparse->indexvals[208] = 74;
  JSparse->indexvals[209] = 75;
  JSparse->indexvals[210] = 75;
  JSparse->indexvals[211] = 76;
  JSparse->indexvals[212] = 76;
  JSparse->indexvals[213] = 77;
  JSparse->indexvals[214] = 77;
  JSparse->indexvals[215] = 78;
  JSparse->indexvals[216] = 78;
  JSparse->indexvals[217] = 79;
  JSparse->indexvals[218] = 79;
  JSparse->indexvals[219] = 80;
  JSparse->indexvals[220] = 72;
  JSparse->indexvals[221] = 80;
  JSparse->indexvals[222] = 81;
  JSparse->indexvals[223] = 82;
  JSparse->indexvals[224] = 82;
  JSparse->indexvals[225] = 83;
  JSparse->indexvals[226] = 83;
  JSparse->indexvals[227] = 84;
  JSparse->indexvals[228] = 84;
  JSparse->indexvals[229] = 85;
  JSparse->indexvals[230] = 85;
  JSparse->indexvals[231] = 86;
  JSparse->indexvals[232] = 86;
  JSparse->indexvals[233] = 87;
  JSparse->indexvals[234] = 87;
  JSparse->indexvals[235] = 88;
  JSparse->indexvals[236] = 88;
  JSparse->indexvals[237] = 89;
  JSparse->indexvals[238] = 81;
  JSparse->indexvals[239] = 89;
  JSparse->indexvals[240] = 90;
  JSparse->indexvals[241] = 91;
  JSparse->indexvals[242] = 91;
  JSparse->indexvals[243] = 92;
  JSparse->indexvals[244] = 92;
  JSparse->indexvals[245] = 93;
  JSparse->indexvals[246] = 93;
  JSparse->indexvals[247] = 94;
  JSparse->indexvals[248] = 94;
  JSparse->indexvals[249] = 95;
  JSparse->indexvals[250] = 95;
  JSparse->indexvals[251] = 96;
  JSparse->indexvals[252] = 96;
  JSparse->indexvals[253] = 97;
  JSparse->indexvals[254] = 97;
  JSparse->indexvals[255] = 98;
  JSparse->indexvals[256] = 90;
  JSparse->indexvals[257] = 98;
  JSparse->indexvals[258] = 99;
  JSparse->indexvals[259] = 100;
  JSparse->indexvals[260] = 100;
  JSparse->indexvals[261] = 101;
  JSparse->indexvals[262] = 101;
  JSparse->indexvals[263] = 102;
  JSparse->indexvals[264] = 102;
  JSparse->indexvals[265] = 103;
  JSparse->indexvals[266] = 103;
  JSparse->indexvals[267] = 104;
  JSparse->indexvals[268] = 104;
  JSparse->indexvals[269] = 105;
  JSparse->indexvals[270] = 105;
  JSparse->indexvals[271] = 106;
  JSparse->indexvals[272] = 106;
  JSparse->indexvals[273] = 107;
  JSparse->indexvals[274] = 99;
  JSparse->indexvals[275] = 107;
  JSparse->indexvals[276] = 108;
  JSparse->indexvals[277] = 109;
  JSparse->indexvals[278] = 109;
  JSparse->indexvals[279] = 110;
  JSparse->indexvals[280] = 110;
  JSparse->indexvals[281] = 111;
  JSparse->indexvals[282] = 111;
  JSparse->indexvals[283] = 112;
  JSparse->indexvals[284] = 112;
  JSparse->indexvals[285] = 113;
  JSparse->indexvals[286] = 113;
  JSparse->indexvals[287] = 114;
  JSparse->indexvals[288] = 114;
  JSparse->indexvals[289] = 115;
  JSparse->indexvals[290] = 115;
  JSparse->indexvals[291] = 116;
  JSparse->indexvals[292] = 108;
  JSparse->indexvals[293] = 116;
  JSparse->indexvals[294] = 117;
  JSparse->indexvals[295] = 118;
  JSparse->indexvals[296] = 118;
  JSparse->indexvals[297] = 119;
  JSparse->indexvals[298] = 119;
  JSparse->indexvals[299] = 120;
  JSparse->indexvals[300] = 120;
  JSparse->indexvals[301] = 121;
  JSparse->indexvals[302] = 121;
  JSparse->indexvals[303] = 122;
  JSparse->indexvals[304] = 122;
  JSparse->indexvals[305] = 123;
  JSparse->indexvals[306] = 123;
  JSparse->indexvals[307] = 124;
  JSparse->indexvals[308] = 124;
  JSparse->indexvals[309] = 125;
  JSparse->indexvals[310] = 117;
  JSparse->indexvals[311] = 125;
  JSparse->indexvals[312] = 126;
  JSparse->indexvals[313] = 127;
  JSparse->indexvals[314] = 127;
  JSparse->indexvals[315] = 128;
  JSparse->indexvals[316] = 128;
  JSparse->indexvals[317] = 129;
  JSparse->indexvals[318] = 129;
  JSparse->indexvals[319] = 130;
  JSparse->indexvals[320] = 130;
  JSparse->indexvals[321] = 131;
  JSparse->indexvals[322] = 131;
  JSparse->indexvals[323] = 132;
  JSparse->indexvals[324] = 132;
  JSparse->indexvals[325] = 133;
  JSparse->indexvals[326] = 133;
  JSparse->indexvals[327] = 134;
  JSparse->indexvals[328] = 126;
  JSparse->indexvals[329] = 134;
  JSparse->indexvals[330] = 135;
  JSparse->indexvals[331] = 136;
  JSparse->indexvals[332] = 136;
  JSparse->indexvals[333] = 137;
  JSparse->indexvals[334] = 137;
  JSparse->indexvals[335] = 138;
  JSparse->indexvals[336] = 138;
  JSparse->indexvals[337] = 139;
  JSparse->indexvals[338] = 139;
  JSparse->indexvals[339] = 140;
  JSparse->indexvals[340] = 140;
  JSparse->indexvals[341] = 141;
  JSparse->indexvals[342] = 141;
  JSparse->indexvals[343] = 142;
  JSparse->indexvals[344] = 142;
  JSparse->indexvals[345] = 143;
  JSparse->indexvals[346] = 135;
  JSparse->indexvals[347] = 143;
  JSparse->indexvals[348] = 144;
  JSparse->indexvals[349] = 145;
  JSparse->indexvals[350] = 145;
  JSparse->indexvals[351] = 146;
  JSparse->indexvals[352] = 146;
  JSparse->indexvals[353] = 147;
  JSparse->indexvals[354] = 147;
  JSparse->indexvals[355] = 148;
  JSparse->indexvals[356] = 148;
  JSparse->indexvals[357] = 149;
  JSparse->indexvals[358] = 149;
  JSparse->indexvals[359] = 150;
  JSparse->indexvals[360] = 150;
  JSparse->indexvals[361] = 151;
  JSparse->indexvals[362] = 151;
  JSparse->indexvals[363] = 152;
  JSparse->indexvals[364] = 144;
  JSparse->indexvals[365] = 152;
  JSparse->indexvals[366] = 153;
  JSparse->indexvals[367] = 154;
  JSparse->indexvals[368] = 154;
  JSparse->indexvals[369] = 155;
  JSparse->indexvals[370] = 155;
  JSparse->indexvals[371] = 156;
  JSparse->indexvals[372] = 156;
  JSparse->indexvals[373] = 157;
  JSparse->indexvals[374] = 157;
  JSparse->indexvals[375] = 158;
  JSparse->indexvals[376] = 158;
  JSparse->indexvals[377] = 159;
  JSparse->indexvals[378] = 159;
  JSparse->indexvals[379] = 160;
  JSparse->indexvals[380] = 160;
  JSparse->indexvals[381] = 161;
  JSparse->indexvals[382] = 153;
  JSparse->indexvals[383] = 161;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 14;
  JSparse->indexptrs[2] = 50;
  JSparse->indexptrs[3] = 54;
  JSparse->indexptrs[4] = 58;
  JSparse->indexptrs[5] = 62;
  JSparse->indexptrs[6] = 66;
  JSparse->indexptrs[7] = 70;
  JSparse->indexptrs[8] = 74;
  JSparse->indexptrs[9] = 78;
  JSparse->indexptrs[10] = 80;
  JSparse->indexptrs[11] = 82;
  JSparse->indexptrs[12] = 84;
  JSparse->indexptrs[13] = 86;
  JSparse->indexptrs[14] = 88;
  JSparse->indexptrs[15] = 90;
  JSparse->indexptrs[16] = 92;
  JSparse->indexptrs[17] = 94;
  JSparse->indexptrs[18] = 96;
  JSparse->indexptrs[19] = 98;
  JSparse->indexptrs[20] = 100;
  JSparse->indexptrs[21] = 102;
  JSparse->indexptrs[22] = 104;
  JSparse->indexptrs[23] = 106;
  JSparse->indexptrs[24] = 108;
  JSparse->indexptrs[25] = 110;
  JSparse->indexptrs[26] = 112;
  JSparse->indexptrs[27] = 114;
  JSparse->indexptrs[28] = 116;
  JSparse->indexptrs[29] = 118;
  JSparse->indexptrs[30] = 120;
  JSparse->indexptrs[31] = 122;
  JSparse->indexptrs[32] = 124;
  JSparse->indexptrs[33] = 126;
  JSparse->indexptrs[34] = 128;
  JSparse->indexptrs[35] = 130;
  JSparse->indexptrs[36] = 132;
  JSparse->indexptrs[37] = 134;
  JSparse->indexptrs[38] = 136;
  JSparse->indexptrs[39] = 138;
  JSparse->indexptrs[40] = 140;
  JSparse->indexptrs[41] = 142;
  JSparse->indexptrs[42] = 144;
  JSparse->indexptrs[43] = 146;
  JSparse->indexptrs[44] = 148;
  JSparse->indexptrs[45] = 150;
  JSparse->indexptrs[46] = 152;
  JSparse->indexptrs[47] = 154;
  JSparse->indexptrs[48] = 156;
  JSparse->indexptrs[49] = 158;
  JSparse->indexptrs[50] = 160;
  JSparse->indexptrs[51] = 162;
  JSparse->indexptrs[52] = 164;
  JSparse->indexptrs[53] = 166;
  JSparse->indexptrs[54] = 168;
  JSparse->indexptrs[55] = 170;
  JSparse->indexptrs[56] = 172;
  JSparse->indexptrs[57] = 174;
  JSparse->indexptrs[58] = 176;
  JSparse->indexptrs[59] = 178;
  JSparse->indexptrs[60] = 180;
  JSparse->indexptrs[61] = 182;
  JSparse->indexptrs[62] = 184;
  JSparse->indexptrs[63] = 186;
  JSparse->indexptrs[64] = 188;
  JSparse->indexptrs[65] = 190;
  JSparse->indexptrs[66] = 192;
  JSparse->indexptrs[67] = 194;
  JSparse->indexptrs[68] = 196;
  JSparse->indexptrs[69] = 198;
  JSparse->indexptrs[70] = 200;
  JSparse->indexptrs[71] = 202;
  JSparse->indexptrs[72] = 204;
  JSparse->indexptrs[73] = 206;
  JSparse->indexptrs[74] = 208;
  JSparse->indexptrs[75] = 210;
  JSparse->indexptrs[76] = 212;
  JSparse->indexptrs[77] = 214;
  JSparse->indexptrs[78] = 216;
  JSparse->indexptrs[79] = 218;
  JSparse->indexptrs[80] = 220;
  JSparse->indexptrs[81] = 222;
  JSparse->indexptrs[82] = 224;
  JSparse->indexptrs[83] = 226;
  JSparse->indexptrs[84] = 228;
  JSparse->indexptrs[85] = 230;
  JSparse->indexptrs[86] = 232;
  JSparse->indexptrs[87] = 234;
  JSparse->indexptrs[88] = 236;
  JSparse->indexptrs[89] = 238;
  JSparse->indexptrs[90] = 240;
  JSparse->indexptrs[91] = 242;
  JSparse->indexptrs[92] = 244;
  JSparse->indexptrs[93] = 246;
  JSparse->indexptrs[94] = 248;
  JSparse->indexptrs[95] = 250;
  JSparse->indexptrs[96] = 252;
  JSparse->indexptrs[97] = 254;
  JSparse->indexptrs[98] = 256;
  JSparse->indexptrs[99] = 258;
  JSparse->indexptrs[100] = 260;
  JSparse->indexptrs[101] = 262;
  JSparse->indexptrs[102] = 264;
  JSparse->indexptrs[103] = 266;
  JSparse->indexptrs[104] = 268;
  JSparse->indexptrs[105] = 270;
  JSparse->indexptrs[106] = 272;
  JSparse->indexptrs[107] = 274;
  JSparse->indexptrs[108] = 276;
  JSparse->indexptrs[109] = 278;
  JSparse->indexptrs[110] = 280;
  JSparse->indexptrs[111] = 282;
  JSparse->indexptrs[112] = 284;
  JSparse->indexptrs[113] = 286;
  JSparse->indexptrs[114] = 288;
  JSparse->indexptrs[115] = 290;
  JSparse->indexptrs[116] = 292;
  JSparse->indexptrs[117] = 294;
  JSparse->indexptrs[118] = 296;
  JSparse->indexptrs[119] = 298;
  JSparse->indexptrs[120] = 300;
  JSparse->indexptrs[121] = 302;
  JSparse->indexptrs[122] = 304;
  JSparse->indexptrs[123] = 306;
  JSparse->indexptrs[124] = 308;
  JSparse->indexptrs[125] = 310;
  JSparse->indexptrs[126] = 312;
  JSparse->indexptrs[127] = 314;
  JSparse->indexptrs[128] = 316;
  JSparse->indexptrs[129] = 318;
  JSparse->indexptrs[130] = 320;
  JSparse->indexptrs[131] = 322;
  JSparse->indexptrs[132] = 324;
  JSparse->indexptrs[133] = 326;
  JSparse->indexptrs[134] = 328;
  JSparse->indexptrs[135] = 330;
  JSparse->indexptrs[136] = 332;
  JSparse->indexptrs[137] = 334;
  JSparse->indexptrs[138] = 336;
  JSparse->indexptrs[139] = 338;
  JSparse->indexptrs[140] = 340;
  JSparse->indexptrs[141] = 342;
  JSparse->indexptrs[142] = 344;
  JSparse->indexptrs[143] = 346;
  JSparse->indexptrs[144] = 348;
  JSparse->indexptrs[145] = 350;
  JSparse->indexptrs[146] = 352;
  JSparse->indexptrs[147] = 354;
  JSparse->indexptrs[148] = 356;
  JSparse->indexptrs[149] = 358;
  JSparse->indexptrs[150] = 360;
  JSparse->indexptrs[151] = 362;
  JSparse->indexptrs[152] = 364;
  JSparse->indexptrs[153] = 366;
  JSparse->indexptrs[154] = 368;
  JSparse->indexptrs[155] = 370;
  JSparse->indexptrs[156] = 372;
  JSparse->indexptrs[157] = 374;
  JSparse->indexptrs[158] = 376;
  JSparse->indexptrs[159] = 378;
  JSparse->indexptrs[160] = 380;
  JSparse->indexptrs[161] = 382;
  JSparse->indexptrs[162] = 384;
  JSparse->data[0] = -k[0]*p[0]*w[0]*w[2];
  JSparse->data[1] = p[0]*w[0];
  JSparse->data[2] = -w[0];
  JSparse->data[3] = w[0];
  JSparse->data[4] = -p[0]*w[5];
  JSparse->data[5] = p[0]*w[5];
  JSparse->data[6] = -p[0]*w[6];
  JSparse->data[7] = p[0]*w[6];
  JSparse->data[8] = -p[0]*w[7];
  JSparse->data[9] = p[0]*w[7];
  JSparse->data[10] = -p[0]*w[8];
  JSparse->data[11] = p[0]*w[8];
  JSparse->data[12] = -p[0]*w[9];
  JSparse->data[13] = p[0]*w[9];
  JSparse->data[14] = p[1]*dwdx[0]*-2.0;
  JSparse->data[15] = p[1]*dwdx[0];
  JSparse->data[16] = p[1]*x[10]*-4.0;
  JSparse->data[17] = p[1]*x[10]*2.0;
  JSparse->data[18] = dwdx[0]*-2.0-p[1]*x[19]*4.0;
  JSparse->data[19] = dwdx[0]+p[1]*x[19]*2.0;
  JSparse->data[20] = p[1]*x[28]*-4.0;
  JSparse->data[21] = p[1]*x[28]*2.0;
  JSparse->data[22] = p[1]*x[37]*-4.0;
  JSparse->data[23] = p[1]*x[37]*2.0;
  JSparse->data[24] = p[1]*x[46]*-4.0;
  JSparse->data[25] = p[1]*x[46]*2.0;
  JSparse->data[26] = p[1]*x[55]*-4.0;
  JSparse->data[27] = p[1]*x[55]*2.0;
  JSparse->data[28] = p[1]*x[64]*-4.0;
  JSparse->data[29] = p[1]*x[64]*2.0;
  JSparse->data[30] = p[1]*x[73]*-4.0;
  JSparse->data[31] = p[1]*x[73]*2.0;
  JSparse->data[32] = p[1]*x[82]*-4.0;
  JSparse->data[33] = p[1]*x[82]*2.0;
  JSparse->data[34] = p[1]*x[91]*-4.0;
  JSparse->data[35] = p[1]*x[91]*2.0;
  JSparse->data[36] = p[1]*x[100]*-4.0;
  JSparse->data[37] = p[1]*x[100]*2.0;
  JSparse->data[38] = p[1]*x[109]*-4.0;
  JSparse->data[39] = p[1]*x[109]*2.0;
  JSparse->data[40] = p[1]*x[118]*-4.0;
  JSparse->data[41] = p[1]*x[118]*2.0;
  JSparse->data[42] = p[1]*x[127]*-4.0;
  JSparse->data[43] = p[1]*x[127]*2.0;
  JSparse->data[44] = p[1]*x[136]*-4.0;
  JSparse->data[45] = p[1]*x[136]*2.0;
  JSparse->data[46] = p[1]*x[145]*-4.0;
  JSparse->data[47] = p[1]*x[145]*2.0;
  JSparse->data[48] = p[1]*x[154]*-4.0;
  JSparse->data[49] = p[1]*x[154]*2.0;
  JSparse->data[50] = -p[2];
  JSparse->data[51] = k[0]*p[2]*w[3];
  JSparse->data[52] = -1.0;
  JSparse->data[53] = k[0]*w[3];
  JSparse->data[54] = -k[1]*p[3]*w[3];
  JSparse->data[55] = p[3]*dwdx[1];
  JSparse->data[56] = -1.0;
  JSparse->data[57] = dwdx[1];
  JSparse->data[58] = -p[3];
  JSparse->data[59] = p[3];
  JSparse->data[60] = -1.0;
  JSparse->data[61] = 1.0;
  JSparse->data[62] = -p[3];
  JSparse->data[63] = p[3];
  JSparse->data[64] = -1.0;
  JSparse->data[65] = 1.0;
  JSparse->data[66] = -p[3];
  JSparse->data[67] = p[3];
  JSparse->data[68] = -1.0;
  JSparse->data[69] = 1.0;
  JSparse->data[70] = -p[3];
  JSparse->data[71] = p[3];
  JSparse->data[72] = -1.0;
  JSparse->data[73] = 1.0;
  JSparse->data[74] = k[1]*p[3]*w[2];
  JSparse->data[75] = -p[3];
  JSparse->data[76] = k[1]*w[2];
  JSparse->data[77] = -1.0;
  JSparse->data[78] = -p[0]*w[0];
  JSparse->data[79] = p[0]*w[0];
  JSparse->data[80] = p[1]*x[1]*-4.0;
  JSparse->data[81] = p[1]*x[1]*2.0;
  JSparse->data[82] = -p[2];
  JSparse->data[83] = k[0]*p[2]*w[3];
  JSparse->data[84] = -p[3];
  JSparse->data[85] = p[3]*2.0;
  JSparse->data[86] = -p[3];
  JSparse->data[87] = p[3];
  JSparse->data[88] = -p[3];
  JSparse->data[89] = p[3];
  JSparse->data[90] = -p[3];
  JSparse->data[91] = p[3];
  JSparse->data[92] = -p[3];
  JSparse->data[93] = p[3];
  JSparse->data[94] = k[1]*p[3]*w[2];
  JSparse->data[95] = -p[3];
  JSparse->data[96] = -p[0]*w[0];
  JSparse->data[97] = p[0]*w[0];
  JSparse->data[98] = p[1]*x[1]*-4.0;
  JSparse->data[99] = p[1]*x[1]*2.0;
  JSparse->data[100] = -p[2];
  JSparse->data[101] = k[0]*p[2]*w[3];
  JSparse->data[102] = -p[3];
  JSparse->data[103] = p[3]*2.0;
  JSparse->data[104] = -p[3];
  JSparse->data[105] = p[3];
  JSparse->data[106] = -p[3];
  JSparse->data[107] = p[3];
  JSparse->data[108] = -p[3];
  JSparse->data[109] = p[3];
  JSparse->data[110] = -p[3];
  JSparse->data[111] = p[3];
  JSparse->data[112] = k[1]*p[3]*w[2];
  JSparse->data[113] = -p[3];
  JSparse->data[114] = -p[0]*w[0];
  JSparse->data[115] = p[0]*w[0];
  JSparse->data[116] = p[1]*x[1]*-4.0;
  JSparse->data[117] = p[1]*x[1]*2.0;
  JSparse->data[118] = -p[2];
  JSparse->data[119] = k[0]*p[2]*w[3];
  JSparse->data[120] = -p[3];
  JSparse->data[121] = p[3]*2.0;
  JSparse->data[122] = -p[3];
  JSparse->data[123] = p[3];
  JSparse->data[124] = -p[3];
  JSparse->data[125] = p[3];
  JSparse->data[126] = -p[3];
  JSparse->data[127] = p[3];
  JSparse->data[128] = -p[3];
  JSparse->data[129] = p[3];
  JSparse->data[130] = k[1]*p[3]*w[2];
  JSparse->data[131] = -p[3];
  JSparse->data[132] = -p[0]*w[0];
  JSparse->data[133] = p[0]*w[0];
  JSparse->data[134] = p[1]*x[1]*-4.0;
  JSparse->data[135] = p[1]*x[1]*2.0;
  JSparse->data[136] = -p[2];
  JSparse->data[137] = k[0]*p[2]*w[3];
  JSparse->data[138] = -p[3];
  JSparse->data[139] = p[3]*2.0;
  JSparse->data[140] = -p[3];
  JSparse->data[141] = p[3];
  JSparse->data[142] = -p[3];
  JSparse->data[143] = p[3];
  JSparse->data[144] = -p[3];
  JSparse->data[145] = p[3];
  JSparse->data[146] = -p[3];
  JSparse->data[147] = p[3];
  JSparse->data[148] = k[1]*p[3]*w[2];
  JSparse->data[149] = -p[3];
  JSparse->data[150] = -p[0]*w[0];
  JSparse->data[151] = p[0]*w[0];
  JSparse->data[152] = p[1]*x[1]*-4.0;
  JSparse->data[153] = p[1]*x[1]*2.0;
  JSparse->data[154] = -p[2];
  JSparse->data[155] = k[0]*p[2]*w[3];
  JSparse->data[156] = -p[3];
  JSparse->data[157] = p[3]*2.0;
  JSparse->data[158] = -p[3];
  JSparse->data[159] = p[3];
  JSparse->data[160] = -p[3];
  JSparse->data[161] = p[3];
  JSparse->data[162] = -p[3];
  JSparse->data[163] = p[3];
  JSparse->data[164] = -p[3];
  JSparse->data[165] = p[3];
  JSparse->data[166] = k[1]*p[3]*w[2];
  JSparse->data[167] = -p[3];
  JSparse->data[168] = -p[0]*w[0];
  JSparse->data[169] = p[0]*w[0];
  JSparse->data[170] = p[1]*x[1]*-4.0;
  JSparse->data[171] = p[1]*x[1]*2.0;
  JSparse->data[172] = -p[2];
  JSparse->data[173] = k[0]*p[2]*w[3];
  JSparse->data[174] = -p[3];
  JSparse->data[175] = p[3]*2.0;
  JSparse->data[176] = -p[3];
  JSparse->data[177] = p[3];
  JSparse->data[178] = -p[3];
  JSparse->data[179] = p[3];
  JSparse->data[180] = -p[3];
  JSparse->data[181] = p[3];
  JSparse->data[182] = -p[3];
  JSparse->data[183] = p[3];
  JSparse->data[184] = k[1]*p[3]*w[2];
  JSparse->data[185] = -p[3];
  JSparse->data[186] = -p[0]*w[0];
  JSparse->data[187] = p[0]*w[0];
  JSparse->data[188] = p[1]*x[1]*-4.0;
  JSparse->data[189] = p[1]*x[1]*2.0;
  JSparse->data[190] = -p[2];
  JSparse->data[191] = k[0]*p[2]*w[3];
  JSparse->data[192] = -p[3];
  JSparse->data[193] = p[3]*2.0;
  JSparse->data[194] = -p[3];
  JSparse->data[195] = p[3];
  JSparse->data[196] = -p[3];
  JSparse->data[197] = p[3];
  JSparse->data[198] = -p[3];
  JSparse->data[199] = p[3];
  JSparse->data[200] = -p[3];
  JSparse->data[201] = p[3];
  JSparse->data[202] = k[1]*p[3]*w[2];
  JSparse->data[203] = -p[3];
  JSparse->data[204] = -p[0]*w[0];
  JSparse->data[205] = p[0]*w[0];
  JSparse->data[206] = p[1]*x[1]*-4.0;
  JSparse->data[207] = p[1]*x[1]*2.0;
  JSparse->data[208] = -p[2];
  JSparse->data[209] = k[0]*p[2]*w[3];
  JSparse->data[210] = -p[3];
  JSparse->data[211] = p[3]*2.0;
  JSparse->data[212] = -p[3];
  JSparse->data[213] = p[3];
  JSparse->data[214] = -p[3];
  JSparse->data[215] = p[3];
  JSparse->data[216] = -p[3];
  JSparse->data[217] = p[3];
  JSparse->data[218] = -p[3];
  JSparse->data[219] = p[3];
  JSparse->data[220] = k[1]*p[3]*w[2];
  JSparse->data[221] = -p[3];
  JSparse->data[222] = -p[0]*w[0];
  JSparse->data[223] = p[0]*w[0];
  JSparse->data[224] = p[1]*x[1]*-4.0;
  JSparse->data[225] = p[1]*x[1]*2.0;
  JSparse->data[226] = -p[2];
  JSparse->data[227] = k[0]*p[2]*w[3];
  JSparse->data[228] = -p[3];
  JSparse->data[229] = p[3]*2.0;
  JSparse->data[230] = -p[3];
  JSparse->data[231] = p[3];
  JSparse->data[232] = -p[3];
  JSparse->data[233] = p[3];
  JSparse->data[234] = -p[3];
  JSparse->data[235] = p[3];
  JSparse->data[236] = -p[3];
  JSparse->data[237] = p[3];
  JSparse->data[238] = k[1]*p[3]*w[2];
  JSparse->data[239] = -p[3];
  JSparse->data[240] = -p[0]*w[0];
  JSparse->data[241] = p[0]*w[0];
  JSparse->data[242] = p[1]*x[1]*-4.0;
  JSparse->data[243] = p[1]*x[1]*2.0;
  JSparse->data[244] = -p[2];
  JSparse->data[245] = k[0]*p[2]*w[3];
  JSparse->data[246] = -p[3];
  JSparse->data[247] = p[3]*2.0;
  JSparse->data[248] = -p[3];
  JSparse->data[249] = p[3];
  JSparse->data[250] = -p[3];
  JSparse->data[251] = p[3];
  JSparse->data[252] = -p[3];
  JSparse->data[253] = p[3];
  JSparse->data[254] = -p[3];
  JSparse->data[255] = p[3];
  JSparse->data[256] = k[1]*p[3]*w[2];
  JSparse->data[257] = -p[3];
  JSparse->data[258] = -p[0]*w[0];
  JSparse->data[259] = p[0]*w[0];
  JSparse->data[260] = p[1]*x[1]*-4.0;
  JSparse->data[261] = p[1]*x[1]*2.0;
  JSparse->data[262] = -p[2];
  JSparse->data[263] = k[0]*p[2]*w[3];
  JSparse->data[264] = -p[3];
  JSparse->data[265] = p[3]*2.0;
  JSparse->data[266] = -p[3];
  JSparse->data[267] = p[3];
  JSparse->data[268] = -p[3];
  JSparse->data[269] = p[3];
  JSparse->data[270] = -p[3];
  JSparse->data[271] = p[3];
  JSparse->data[272] = -p[3];
  JSparse->data[273] = p[3];
  JSparse->data[274] = k[1]*p[3]*w[2];
  JSparse->data[275] = -p[3];
  JSparse->data[276] = -p[0]*w[0];
  JSparse->data[277] = p[0]*w[0];
  JSparse->data[278] = p[1]*x[1]*-4.0;
  JSparse->data[279] = p[1]*x[1]*2.0;
  JSparse->data[280] = -p[2];
  JSparse->data[281] = k[0]*p[2]*w[3];
  JSparse->data[282] = -p[3];
  JSparse->data[283] = p[3]*2.0;
  JSparse->data[284] = -p[3];
  JSparse->data[285] = p[3];
  JSparse->data[286] = -p[3];
  JSparse->data[287] = p[3];
  JSparse->data[288] = -p[3];
  JSparse->data[289] = p[3];
  JSparse->data[290] = -p[3];
  JSparse->data[291] = p[3];
  JSparse->data[292] = k[1]*p[3]*w[2];
  JSparse->data[293] = -p[3];
  JSparse->data[294] = -p[0]*w[0];
  JSparse->data[295] = p[0]*w[0];
  JSparse->data[296] = p[1]*x[1]*-4.0;
  JSparse->data[297] = p[1]*x[1]*2.0;
  JSparse->data[298] = -p[2];
  JSparse->data[299] = k[0]*p[2]*w[3];
  JSparse->data[300] = -p[3];
  JSparse->data[301] = p[3]*2.0;
  JSparse->data[302] = -p[3];
  JSparse->data[303] = p[3];
  JSparse->data[304] = -p[3];
  JSparse->data[305] = p[3];
  JSparse->data[306] = -p[3];
  JSparse->data[307] = p[3];
  JSparse->data[308] = -p[3];
  JSparse->data[309] = p[3];
  JSparse->data[310] = k[1]*p[3]*w[2];
  JSparse->data[311] = -p[3];
  JSparse->data[312] = -p[0]*w[0];
  JSparse->data[313] = p[0]*w[0];
  JSparse->data[314] = p[1]*x[1]*-4.0;
  JSparse->data[315] = p[1]*x[1]*2.0;
  JSparse->data[316] = -p[2];
  JSparse->data[317] = k[0]*p[2]*w[3];
  JSparse->data[318] = -p[3];
  JSparse->data[319] = p[3]*2.0;
  JSparse->data[320] = -p[3];
  JSparse->data[321] = p[3];
  JSparse->data[322] = -p[3];
  JSparse->data[323] = p[3];
  JSparse->data[324] = -p[3];
  JSparse->data[325] = p[3];
  JSparse->data[326] = -p[3];
  JSparse->data[327] = p[3];
  JSparse->data[328] = k[1]*p[3]*w[2];
  JSparse->data[329] = -p[3];
  JSparse->data[330] = -p[0]*w[0];
  JSparse->data[331] = p[0]*w[0];
  JSparse->data[332] = p[1]*x[1]*-4.0;
  JSparse->data[333] = p[1]*x[1]*2.0;
  JSparse->data[334] = -p[2];
  JSparse->data[335] = k[0]*p[2]*w[3];
  JSparse->data[336] = -p[3];
  JSparse->data[337] = p[3]*2.0;
  JSparse->data[338] = -p[3];
  JSparse->data[339] = p[3];
  JSparse->data[340] = -p[3];
  JSparse->data[341] = p[3];
  JSparse->data[342] = -p[3];
  JSparse->data[343] = p[3];
  JSparse->data[344] = -p[3];
  JSparse->data[345] = p[3];
  JSparse->data[346] = k[1]*p[3]*w[2];
  JSparse->data[347] = -p[3];
  JSparse->data[348] = -p[0]*w[0];
  JSparse->data[349] = p[0]*w[0];
  JSparse->data[350] = p[1]*x[1]*-4.0;
  JSparse->data[351] = p[1]*x[1]*2.0;
  JSparse->data[352] = -p[2];
  JSparse->data[353] = k[0]*p[2]*w[3];
  JSparse->data[354] = -p[3];
  JSparse->data[355] = p[3]*2.0;
  JSparse->data[356] = -p[3];
  JSparse->data[357] = p[3];
  JSparse->data[358] = -p[3];
  JSparse->data[359] = p[3];
  JSparse->data[360] = -p[3];
  JSparse->data[361] = p[3];
  JSparse->data[362] = -p[3];
  JSparse->data[363] = p[3];
  JSparse->data[364] = k[1]*p[3]*w[2];
  JSparse->data[365] = -p[3];
  JSparse->data[366] = -p[0]*w[0];
  JSparse->data[367] = p[0]*w[0];
  JSparse->data[368] = p[1]*x[1]*-4.0;
  JSparse->data[369] = p[1]*x[1]*2.0;
  JSparse->data[370] = -p[2];
  JSparse->data[371] = k[0]*p[2]*w[3];
  JSparse->data[372] = -p[3];
  JSparse->data[373] = p[3]*2.0;
  JSparse->data[374] = -p[3];
  JSparse->data[375] = p[3];
  JSparse->data[376] = -p[3];
  JSparse->data[377] = p[3];
  JSparse->data[378] = -p[3];
  JSparse->data[379] = p[3];
  JSparse->data[380] = -p[3];
  JSparse->data[381] = p[3];
  JSparse->data[382] = k[1]*p[3]*w[2];
  JSparse->data[383] = -p[3];
}

} // namespace model_model_jakstat_adjoint_o2

} // namespace amici

