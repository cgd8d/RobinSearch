#ifndef QueueContainers_hpp
#define QueueContainers_hpp

/*
The purpose of this file is to hold queued
data precomputed for use in adding new primes.
This is in two parts:
* PrimeQueueContainer, to hold primes that have
not been added yet.
* TmpProdContainer, to hold products of bunches
of values.

These containers are special-purpose because
they are essential for performance and have
access patterns which I may be able to take
advantage of.  Putting it into a structure
makes it easier to experiment with that.
*/







#endif
