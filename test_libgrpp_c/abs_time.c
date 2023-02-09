/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include "abs_time.h"

#include <stdlib.h>
#include <sys/time.h>

double abs_time()
{
    struct timeval cur_time;
    gettimeofday(&cur_time, NULL);
    return (cur_time.tv_sec * 1000000u + cur_time.tv_usec) / 1.e6;
}
