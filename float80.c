#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// ===== جدول مساعد لحساب طول البتات =====
const uint8_t len8tab[256] = {
    0x00, 0x01, 0x02, 0x02, 0x03, 0x03, 0x03, 0x03, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04,
    0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05, 0x05,
    0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06,
    0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06,
    0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
    0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
    0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
    0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08,
    0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08
};

int Len64(uint64_t x)
{
    int n = 0;
    if (x >= (1ULL << 32)) { x >>= 32; n = 32; }
    if (x >= (1ULL << 16)) { x >>= 16; n += 16; }
    if (x >= (1ULL << 8))  { x >>= 8;  n += 8; }
    return n + len8tab[(uint8_t)x];
}

int LeadingZeros64(uint64_t x) { 
    return x == 0 ? 64 : 64 - Len64(x); 
}

int LeadingZeros128(uint64_t hi, uint64_t lo)
{
    if (hi != 0) return LeadingZeros64(hi);
    return 64 + LeadingZeros64(lo);
}


// ===== De Bruijn Sequence لحساب Trailing Zeros =====
static const uint64_t deBruijn64 = 0x03f79d71b4ca8b09ULL;

static const uint8_t deBruijn64tab[64] = {
     0,  1, 56,  2, 57, 49, 28,  3, 61, 58, 42, 50, 38, 29, 17,  4,
    62, 47, 59, 36, 45, 43, 51, 22, 53, 39, 33, 30, 24, 18, 12,  5,
    63, 55, 48, 27, 60, 41, 37, 16, 46, 35, 44, 21, 52, 32, 23, 11,
    54, 26, 40, 15, 34, 20, 31, 10, 25, 14, 19,  9, 13,  8,  7,  6,
};

int TrailingZeros64(uint64_t x)
{
    if (x == 0) return 64;
    return (int)deBruijn64tab[((x & -x) * deBruijn64) >> (64 - 6)];
}

void Add64(uint64_t x, uint64_t y, uint64_t carry, uint64_t* sum, uint64_t* carryOut)
{
    *sum = x + y + carry;
    *carryOut = ((x & y) | ((x | y) & ~(*sum))) >> 63;
}

void Sub64(uint64_t x, uint64_t y, uint64_t borrow, uint64_t* diff, uint64_t* borrowOut)
{
    *diff = x - y - borrow;
    *borrowOut = ((~x & y) | (~(x ^ y) & *diff)) >> 63;
}

void Mul64(uint64_t x, uint64_t y, uint64_t* hi, uint64_t* lo)
{
    const uint64_t mask32 = (1ULL << 32) - 1;
    uint64_t x0 = x & mask32, x1 = x >> 32;
    uint64_t y0 = y & mask32, y1 = y >> 32;
    uint64_t w0 = x0 * y0;
    uint64_t t = x1 * y0 + (w0 >> 32);
    uint64_t w1 = t & mask32;
    uint64_t w2 = t >> 32;
    w1 += x0 * y1;
    *hi = x1 * y1 + w2 + (w1 >> 32);
    *lo = x * y;
}

void Div64(uint64_t hi, uint64_t lo, uint64_t y, uint64_t* quo, uint64_t* rem)
{
    if (y == 0) { *quo = 0; *rem = 0; return; }
    if (hi == 0) { *quo = lo / y; *rem = lo % y; return; }
    if (y <= hi) { *quo = 0; *rem = 0; return; }

    __uint128_t dividend = ((__uint128_t)hi << 64) | lo;
    *quo = (uint64_t)(dividend / y);
    *rem = (uint64_t)(dividend % y);
}
// ===== Float80 Structure =====
typedef struct {
    uint64_t significand;
    uint16_t exponent;
    uint16_t sign;
} float80_t;

#define FLOAT80_BIAS 16383
#define FLOAT80_EXP_MAX 0x7FFE
#define FLOAT80_EXP_MIN 0x0001
#define FLOAT80_EXP_SPECIAL 0x7FFF
#define FLOAT80_EXP_DENORMAL 0x0000
#define INTEGER_BIT (1ULL << 63)

typedef enum {
    ROUND_TO_NEAREST = 0,
    ROUND_UP = 1,
    ROUND_DOWN = 2,
    ROUND_CHOP = 3
} rounding_mode_t;

__thread rounding_mode_t current_rounding_mode = ROUND_TO_NEAREST;

// ===== القيم الخاصة =====
float80_t float80_zero(int sign)
{
    float80_t result;
    result.significand = 0;
    result.exponent = FLOAT80_EXP_DENORMAL;
    result.sign = sign & 1;
    return result;
}

float80_t float80_nan(void)
{
    float80_t result;
    result.significand = 0xC000000000000001ULL;
    result.exponent = FLOAT80_EXP_SPECIAL;
    result.sign = 0;
    return result;
}

float80_t float80_infinity(int sign)
{
    float80_t result;
    result.significand = INTEGER_BIT;
    result.exponent = FLOAT80_EXP_SPECIAL;
    result.sign = sign & 1;
    return result;
}

// ===== دوال الفحص =====
bool float80_is_zero(float80_t x)
{
    return x.exponent == FLOAT80_EXP_DENORMAL && x.significand == 0;
}

bool float80_is_nan(float80_t x)
{
    return x.exponent == FLOAT80_EXP_SPECIAL && 
           (x.significand & 0x7FFFFFFFFFFFFFFFULL) != 0;
}

bool float80_is_infinity(float80_t x)
{
    return x.exponent == FLOAT80_EXP_SPECIAL && 
           (x.significand & 0x7FFFFFFFFFFFFFFFULL) == 0;
}

bool float80_is_denormal(float80_t x)
{
    return x.exponent == FLOAT80_EXP_DENORMAL && x.significand != 0;
}

bool float80_is_supported(float80_t x)
{
    if (x.exponent == FLOAT80_EXP_DENORMAL) {
        return (x.significand >> 63) == 0;
    }
    return (x.significand >> 63) == 1;
}

// ===== إزاحة 128-bit لليمين مع التقريب =====
void shift_right_128_with_rounding(
    uint64_t* hi, 
    uint64_t* lo, 
    int shift, 
    int sign)
{
    if (shift == 0) return;
    
    if (shift >= 128) {
        bool round_away = (current_rounding_mode == ROUND_UP && !sign) || 
                         (current_rounding_mode == ROUND_DOWN && sign);
        if (round_away && (*hi != 0 || *lo != 0)) {
            *hi = 0;
            *lo = 1;
        } else {
            *hi = 0;
            *lo = 0;
        }
        return;
    }

    int guard_bit = 0;
    uint64_t rest_hi = 0, rest_lo = 0;

    if (shift <= 64) {
        if (shift > 0) {
            guard_bit = (*lo >> (shift - 1)) & 1;
            if (shift > 1) {
                rest_lo = *lo & ((1ULL << (shift - 1)) - 1);
            }
        }
        
        if (shift == 64) {
            *lo = *hi;
            *hi = 0;
        } else {
            uint64_t new_lo = (*hi << (64 - shift)) | (*lo >> shift);
            uint64_t new_hi = *hi >> shift;
            *hi = new_hi;
            *lo = new_lo;
        }
    } else {
        int shift2 = shift - 64;
        if (shift2 > 0) {
            guard_bit = (*hi >> (shift2 - 1)) & 1;
            if (shift2 > 1) {
                rest_hi = *hi & ((1ULL << (shift2 - 1)) - 1);
            }
        }
        rest_lo = *lo;
        
        *lo = *hi >> shift2;
        *hi = 0;
    }

    if (guard_bit == 0 && rest_hi == 0 && rest_lo == 0) {
        return;
    }

    bool round_away = (current_rounding_mode == ROUND_UP && !sign) || 
                     (current_rounding_mode == ROUND_DOWN && sign);
    
    if (round_away) {
        uint64_t carry;
        Add64(*lo, 1, 0, lo, &carry);
        if (carry) (*hi)++;
    } else if (current_rounding_mode == ROUND_TO_NEAREST && guard_bit) {
        if (rest_hi != 0 || rest_lo != 0) {
            uint64_t carry;
            Add64(*lo, 1, 0, lo, &carry);
            if (carry) (*hi)++;
        } else if (*lo & 1) {
            uint64_t carry;
            Add64(*lo, 1, 0, lo, &carry);
            if (carry) (*hi)++;
        }
    }
}

// ===== التطبيع =====
float80_t normalize_to_float80(uint64_t signif_hi, uint64_t signif_lo, int true_exponent, int sign)
{
    float80_t result;
    result.sign = sign;
    
    // التعامل مع الصفر
    if (signif_hi == 0 && signif_lo == 0) {
        result.significand = 0;
        result.exponent = FLOAT80_EXP_DENORMAL;
        return result;
    }
    
    // حساب عدد الأصفار البادئة
    int leading_zeros = LeadingZeros128(signif_hi, signif_lo);
    
    // إزاحة لليسار لجعل bit 127 = 1
    if (leading_zeros > 0) {
        if (leading_zeros < 64) {
            signif_hi = (signif_hi << leading_zeros) | (signif_lo >> (64 - leading_zeros));
            signif_lo = signif_lo << leading_zeros;
        } else {
            signif_hi = signif_lo << (leading_zeros - 64);
            signif_lo = 0;
        }
        true_exponent -= leading_zeros;
    }
    
    // التحقق من Overflow
    if (true_exponent > (int)(FLOAT80_EXP_MAX - FLOAT80_BIAS)) {
        if (signif_hi == INTEGER_BIT && signif_lo == 0) {
            return float80_infinity(sign);
        } else if ((current_rounding_mode == ROUND_UP && sign) || 
                   (current_rounding_mode == ROUND_DOWN && !sign) || 
                   current_rounding_mode == ROUND_CHOP) {
            result.exponent = FLOAT80_EXP_MAX;
            result.significand = 0xFFFFFFFFFFFFFFFFULL;
            return result;
        } else {
            return float80_infinity(sign);
        }
    }
    
    // التحقق من Underflow
    if (true_exponent < (int)(FLOAT80_EXP_MIN - FLOAT80_BIAS)) {
        int denormal_shift = (int)(FLOAT80_EXP_MIN - FLOAT80_BIAS) - true_exponent;
        shift_right_128_with_rounding(&signif_hi, &signif_lo, denormal_shift, sign);
        result.exponent = FLOAT80_EXP_DENORMAL;
        
        // بعد الإزاحة، نأخذ الـ 64-bit العليا
        shift_right_128_with_rounding(&signif_hi, &signif_lo, 64, sign);
        result.significand = signif_lo;
        
        // التأكد من أن bit 63 = 0 للـ denormal
        result.significand &= ~INTEGER_BIT;
        return result;
    }
    
    // الحالة العادية: تقريب إلى 64-bit
    shift_right_128_with_rounding(&signif_hi, &signif_lo, 64, sign);
    
    // معالجة حالة تجاوز التقريب (حصلنا على 65-bit)
    if (signif_hi != 0) {
        signif_lo = (signif_hi << 63) | (signif_lo >> 1);
        true_exponent++;
        
        // إعادة التحقق من Overflow
        if (true_exponent > (int)(FLOAT80_EXP_MAX - FLOAT80_BIAS)) {
            return float80_infinity(sign);
        }
    }
    
    result.significand = signif_lo;
    result.exponent = true_exponent + FLOAT80_BIAS;
    
    // التأكد من أن bit 63 = 1
    if (result.exponent != FLOAT80_EXP_DENORMAL) {
        result.significand |= INTEGER_BIT;
    }
    
    return result;
}

// ===== التحويل من/إلى int =====
float80_t float80_from_int(int64_t value)
{
    if (value == 0) {
        return float80_zero(0);
    }
    
    int sign = 0;
    uint64_t abs_value = value;
    
    if (value < 0) {
        sign = 1;
        abs_value = -(uint64_t)value;
    }
    
    // الرقم يبدأ كـ 1.xxxxx × 2^63
    // نحتاج لتطبيعه
    int leading_zeros = LeadingZeros64(abs_value);
    abs_value <<= leading_zeros;
    int true_exponent = 63 - leading_zeros;
    
    float80_t result;
    result.significand = abs_value;
    result.exponent = true_exponent + FLOAT80_BIAS;
    result.sign = sign;
    
    return result;
}

int64_t float80_to_int(float80_t x)
{
    if (!float80_is_supported(x)) return INT64_MIN;
    if (float80_is_nan(x)) return INT64_MIN;
    if (float80_is_infinity(x)) return INT64_MIN;
    
    int true_exponent = (int)x.exponent - FLOAT80_BIAS;
    
    // الأس أكبر من 63 = overflow
    if (true_exponent > 63) return INT64_MIN;
    
    // الأس سالب = النتيجة بين 0 و 1
    if (true_exponent < 0) return 0;
    
    uint64_t result = x.significand;
    
    if (true_exponent < 63) {
        result >>= (63 - true_exponent);
    }
    
    if (x.sign) {
        return result == (UINT64_C(1) << 63) ? INT64_MIN : -(int64_t)result;
    }
    
    return (int64_t)result;
}

// ===== التحويل من/إلى double =====
float80_t float80_from_double(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(double));
    
    uint64_t sign = (bits >> 63) & 1;
    uint64_t exp64 = (bits >> 52) & 0x7FF;
    uint64_t signif64 = bits & 0xFFFFFFFFFFFFFULL;
    
    float80_t result;
    result.sign = sign;
    
    // Zero
    if (exp64 == 0 && signif64 == 0) {
        result.significand = 0;
        result.exponent = FLOAT80_EXP_DENORMAL;
        return result;
    }
    
    // Infinity أو NaN
    if (exp64 == 0x7FF) {
        result.exponent = FLOAT80_EXP_SPECIAL;
        if (signif64 == 0) {
            result.significand = INTEGER_BIT;  // Infinity
        } else {
            result.significand = 0xC000000000000001ULL;  // NaN
        }
        return result;
    }
    
    // Denormal في double
    if (exp64 == 0) {
        // تحويل denormal إلى normalized في float80
        uint64_t signif80 = signif64 << 11;
        int leading_zeros = LeadingZeros64(signif80);
        signif80 <<= leading_zeros;
        int true_exponent = (1 - 1023) - leading_zeros;
        
        result.significand = signif80;
        result.exponent = true_exponent + FLOAT80_BIAS;
        return result;
    }
    
    // Normal
    int true_exponent = (int)exp64 - 1023;
    
    result.significand = (signif64 << 11) | INTEGER_BIT;
    result.exponent = true_exponent + FLOAT80_BIAS;
    
    return result;
}

double float80_to_double(float80_t x)
{
    if (!float80_is_supported(x)) return NAN;
    
    // Zero
    if (float80_is_zero(x)) {
        return x.sign ? -0.0 : 0.0;
    }
    
    // Infinity
    if (float80_is_infinity(x)) {
        return x.sign ? -INFINITY : INFINITY;
    }
    
    // NaN
    if (float80_is_nan(x)) {
        return NAN;
    }
    
    int true_exponent = (int)x.exponent - FLOAT80_BIAS;
    if (x.exponent == FLOAT80_EXP_DENORMAL) {
        true_exponent = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    
    int exp64 = true_exponent + 1023;
    
    // Overflow
    if (exp64 >= 0x7FF) {
        return x.sign ? -INFINITY : INFINITY;
    }
    
    // Underflow إلى denormal في double
    uint64_t signif64;
    if (exp64 <= 0) {
        int shift = 11 - exp64;
        uint64_t hi = 0, lo = x.significand;
        shift_right_128_with_rounding(&hi, &lo, shift, x.sign);
        signif64 = lo & 0xFFFFFFFFFFFFFULL;
        exp64 = 0;
    } else {
        uint64_t hi = 0, lo = x.significand;
        shift_right_128_with_rounding(&hi, &lo, 11, x.sign);
        
        // معالجة حالة تجاوز التقريب
        if (lo & (1ULL << 53)) {
            lo >>= 1;
            exp64++;
            if (exp64 >= 0x7FF) {
                return x.sign ? -INFINITY : INFINITY;
            }
        }
        
        signif64 = lo & 0xFFFFFFFFFFFFFULL;
    }
    
    uint64_t bits = ((uint64_t)x.sign << 63) | 
                    ((uint64_t)exp64 << 52) | 
                    signif64;
    
    double result;
    memcpy(&result, &bits, sizeof(double));
    return result;
}

// ===== العمليات الحسابية =====
float80_t float80_neg(float80_t x)
{
    x.sign = !x.sign;
    return x;
}

float80_t float80_abs(float80_t x)
{
    x.sign = 0;
    return x;
}

float80_t float80_add(float80_t a, float80_t b)
{
    if (!float80_is_supported(a) || !float80_is_supported(b)) {
        return float80_nan();
    }
    if (float80_is_nan(a)) return a;
    if (float80_is_nan(b)) return b;
    
    // ترتيب حسب الأس الأكبر
    if (a.exponent < b.exponent) {
        float80_t temp = a;
        a = b;
        b = temp;
    }
    
    bool flipped = false;
    if (a.sign) {
        a.sign = !a.sign;
        b.sign = !b.sign;
        flipped = true;
    }
    
    uint64_t a_hi = 0, a_lo = a.significand;
    uint64_t b_hi = 0, b_lo = b.significand;
    
    // محاذاة الأس
    int exp_diff = a.exponent - b.exponent;
    shift_right_128_with_rounding(&b_hi, &b_lo, exp_diff, b.sign ^ flipped);
    
    int true_exp_a = (int)a.exponent - FLOAT80_BIAS;
    if (a.exponent == FLOAT80_EXP_DENORMAL) {
        true_exp_a = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    
    int sign = a.sign;
    uint64_t result_hi, result_lo;
    
    if (!b.sign) {
        // جمع
        uint64_t carry1, carry2;
        Add64(a_lo, b_lo, 0, &result_lo, &carry1);
        Add64(a_hi, b_hi, carry1, &result_hi, &carry2);
        
        if (!float80_is_infinity(a) && carry1) {
            shift_right_128_with_rounding(&result_hi, &result_lo, 1, sign);
            true_exp_a++;
        }
    } else {
        // طرح
        if (float80_is_infinity(a) && float80_is_infinity(b)) {
            return float80_nan();
        }
        
        if (current_rounding_mode == ROUND_CHOP && 
            b_lo == 0 && b_hi == 0 && b.significand != 0) {
            b_lo = 1;
        }
        
        if (float80_is_infinity(a)) {
            b_hi = 0;
            b_lo = 0;
        }
        
        if (a_hi > b_hi || (a_hi == b_hi && a_lo >= b_lo)) {
            uint64_t borrow1, borrow2;
            Sub64(a_lo, b_lo, 0, &result_lo, &borrow1);
            Sub64(a_hi, b_hi, borrow1, &result_hi, &borrow2);
        } else {
            uint64_t borrow1, borrow2;
            Sub64(b_lo, a_lo, 0, &result_lo, &borrow1);
            Sub64(b_hi, a_hi, borrow1, &result_hi, &borrow2);
            sign = 1;
        }
        
        if (result_hi == 0 && result_lo == 0) {
            if (a_lo != 0 && current_rounding_mode == ROUND_DOWN) {
                return float80_zero(1);
            }
            return float80_zero(0);
        }
    }
    
    if (flipped) sign = !sign;
    
    // الجمع/الطرح ينتج 64-bit (أو 65-bit مع carry).
    // الدالة normalize_to_float80 تتعامل مع 128-bit (bit 127 = 1).
    // للنتائج التي تناسب 64-bit (result_hi == 0)، نضيف 64 للأس للتعويض.
    if (result_hi == 0 && result_lo != 0) {
        true_exp_a += 64;
    }
    
    float80_t result = normalize_to_float80(result_hi, result_lo, true_exp_a, sign);
    return result;
}

float80_t float80_sub(float80_t a, float80_t b)
{
    return float80_add(a, float80_neg(b));
}

float80_t float80_mul(float80_t a, float80_t b)
{
    if (!float80_is_supported(a) || !float80_is_supported(b)) {
        return float80_nan();
    }
    if (float80_is_nan(a)) return a;
    if (float80_is_nan(b)) return b;
    
    if (float80_is_infinity(a) || float80_is_infinity(b)) {
        if (float80_is_zero(a) || float80_is_zero(b)) {
            return float80_nan();
        }
        return float80_infinity(a.sign ^ b.sign);
    }
    
    int true_exp_a = (int)a.exponent - FLOAT80_BIAS;
    int true_exp_b = (int)b.exponent - FLOAT80_BIAS;
    
    if (a.exponent == FLOAT80_EXP_DENORMAL) {
        true_exp_a = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    if (b.exponent == FLOAT80_EXP_DENORMAL) {
        true_exp_b = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    
    int result_exp = true_exp_a + true_exp_b + 1;
    
    uint64_t result_hi, result_lo;
    Mul64(a.significand, b.significand, &result_hi, &result_lo);
    
    float80_t result = normalize_to_float80(result_hi, result_lo, 
                                           result_exp, a.sign ^ b.sign);
    return result;
}

float80_t float80_div(float80_t a, float80_t b)
{
    if (!float80_is_supported(a) || !float80_is_supported(b)) {
        return float80_nan();
    }
    if (float80_is_nan(a)) return a;
    if (float80_is_nan(b)) return b;
    
    if (float80_is_infinity(a)) {
        if (float80_is_infinity(b)) return float80_nan();
        return float80_infinity(a.sign ^ b.sign);
    }
    
    if (float80_is_infinity(b)) {
        return float80_zero(a.sign ^ b.sign);
    }
    
    if (float80_is_zero(b)) {
        if (float80_is_zero(a)) return float80_nan();
        return float80_infinity(a.sign ^ b.sign);
    }
    
    if (float80_is_zero(a)) {
        return float80_zero(a.sign ^ b.sign);
    }
    
    // القسمة: (a.significand × 2^63) / b.significand
    // نوسع المقسوم بـ 64 بت إضافية للحصول على دقة عالية
    uint64_t dividend_hi = a.significand >> 1;
    uint64_t dividend_lo = a.significand << 63;
    uint64_t b_signif = b.significand;
    
    uint64_t quotient, remainder;
    Div64(dividend_hi, dividend_lo, b_signif, &quotient, &remainder);
    
    // نستخدم الباقي لزيادة الدقة: نضرب الباقي في 2^64 ونقسم مرة أخرى
    uint64_t signif_hi = quotient;
    uint64_t signif_lo = 0;
    
    if (remainder != 0) {
        uint64_t extra_quotient, extra_remainder;
        Div64(remainder, 0, b_signif, &extra_quotient, &extra_remainder);
        signif_lo = extra_quotient;
    }
    
    int true_exp_a = (int)a.exponent - FLOAT80_BIAS;
    int true_exp_b = (int)b.exponent - FLOAT80_BIAS;
    
    if (a.exponent == FLOAT80_EXP_DENORMAL) {
        true_exp_a = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    if (b.exponent == FLOAT80_EXP_DENORMAL) {
        true_exp_b = (int)FLOAT80_EXP_MIN - FLOAT80_BIAS;
    }
    
    int result_exp = true_exp_a - true_exp_b;
    
    float80_t result = normalize_to_float80(signif_hi, signif_lo, 
                                           result_exp, a.sign ^ b.sign);
    return result;
}

// ===== المقارنة =====
bool float80_equals(float80_t a, float80_t b)
{
    if (!float80_is_supported(a) || !float80_is_supported(b)) return false;
    if (float80_is_nan(a) || float80_is_nan(b)) return false;
    
    if (float80_is_zero(a) && float80_is_zero(b)) return true;
    
    return a.sign == b.sign && 
           a.exponent == b.exponent && 
           a.significand == b.significand;
}

bool float80_less_than(float80_t a, float80_t b)
{
    if (!float80_is_supported(a) || !float80_is_supported(b)) return false;
    if (float80_is_nan(a) || float80_is_nan(b)) return false;
    
    if (float80_is_infinity(a) && float80_is_infinity(b) && a.sign == b.sign) {
        return false;
    }
    
    if (float80_is_zero(a) && float80_is_zero(b)) return false;
    
    float80_t diff = float80_sub(a, b);
    return diff.sign == 1 && !float80_is_zero(diff);
}

bool float80_greater_than(float80_t a, float80_t b)
{
    return !float80_less_than(a, b) && !float80_equals(a, b);
}

// ===== الدوال الرياضية =====
float80_t float80_sqrt(float80_t x)
{
    if (float80_is_zero(x)) return x;
    if (float80_is_nan(x) || x.sign) return float80_nan();
    
    float80_t guess = x;
    int true_exp = (int)guess.exponent - FLOAT80_BIAS;
    guess.exponent = FLOAT80_BIAS + (true_exp / 2);
    
    float80_t two = float80_from_int(2);
    float80_t old_guess;
    int iterations = 0;
    
    do {
        old_guess = guess;
        guess = float80_div(
            float80_add(guess, float80_div(x, guess)), 
            two
        );
    } while (!float80_equals(guess, old_guess) && iterations++ < 100);
    
    return guess;
}

// ===== الطباعة =====
void float80_print(float80_t x)
{
    printf("%+.19e\n", float80_to_double(x));
}

void float80_print_bits(float80_t x)
{
    printf("sign=%u, exponent=0x%04X (%d), significand=0x%016llX\n",
           x.sign, x.exponent, (int)x.exponent - FLOAT80_BIAS,
           (unsigned long long)x.significand);
}

// ===== الاختبار =====
int main(void)
{
    printf("=== Float80 Library (64-bit only, Fixed Version) ===\n\n");
    
    float80_t a = float80_from_double(10.5);
    float80_t b = float80_from_double(3.25);
    
    printf("a = "); float80_print(a);
    float80_print_bits(a);
    
    printf("\nb = "); float80_print(b);
    float80_print_bits(b);
    
    printf("\n--- Addition ---\n");
    float80_t sum = float80_add(a, b);
    printf("a + b = "); float80_print(sum);
    float80_print_bits(sum);
    printf("Expected: 13.75\n");
    
    printf("\n--- Subtraction ---\n");
    float80_t diff = float80_sub(a, b);
    printf("a - b = "); float80_print(diff);
    float80_print_bits(diff);
    printf("Expected: 7.25\n");
    
    printf("\n--- Multiplication ---\n");
    float80_t prod = float80_mul(a, b);
    printf("a * b = "); float80_print(prod);
    float80_print_bits(prod);
    printf("Expected: 34.125\n");
    
    printf("\n--- Division ---\n");
    float80_t quot = float80_div(a, b);
    printf("a / b = "); float80_print(quot);
    float80_print_bits(quot);
    printf("Expected: 3.230769...\n");
    
    printf("\n--- Square Root ---\n");
    float80_t sq = float80_sqrt(a);
    printf("sqrt(a) = "); float80_print(sq);
    float80_print_bits(sq);
    printf("Expected: 3.240370...\n");
    
    printf("\n--- Comparison ---\n");
    printf("a == b: %d\n", float80_equals(a, b));
    printf("a > b: %d\n", float80_greater_than(a, b));
    printf("a < b: %d\n", float80_less_than(a, b));
    
    return 0;
}
