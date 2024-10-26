//
// Created by Sam Thompson on 14/03/2019.
//

#ifndef RNG_TESTING_RANDOM_CONVERT_H
#define RNG_TESTING_RANDOM_CONVERT_H

/**
 * @brief Rotate the value left a specified number of times.
 * @param value the value to rotate
 * @param bits the number of bits to shift
 * @return the rotated and shifted value
 */
static inline uint64_t rotl(const uint64_t &value, uint8_t bits)
{
    return (value << bits) | (value >> (64 - bits));
}

/**
 * @brief Rotate the value right a specified number of times.
 * @param value the value to rotate
 * @param bits the number of bits to shift
 * @return the rotated and shifted value
 */
static inline uint64_t rotr(const uint64_t &value, uint8_t rot)
{
    constexpr uint8_t bits = sizeof(uint64_t) * 8;
    constexpr uint8_t mask = bits - 1;
    return (value >> rot) | (value << ((- rot) & mask));
}

/**
 * @brief Converts a uniform integer on [0, UINT64_MAX) to a double
 * @param i the integer to convert
 * @return double equivalent of the integer
 */
static inline double intToDouble(const uint64_t &i)
{
    const union
    {
        uint64_t i;
        double d;
    } u = {.i = UINT64_C(0x3FF) << 52 | i >> 12};
    return u.d - 1.0;
}
#endif //RNG_TESTING_RANDOM_CONVERT_H
