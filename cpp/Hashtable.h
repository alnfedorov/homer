
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@gmail.com>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.


#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <limits.h>


#ifndef HASHTABLE_H
#define HASHTABLE_H


#define EMPTY_INT -1080706050
#define EMPTY_DOUBLE -1.23e100
#define EMPTY_DOUBLE_CHECK -1.20e100

#define DEFAULT_SIZE 100000


class LinkedListItem;
class LinkedList {
public:
	LinkedListItem* firstItem;
	LinkedListItem* lastItem;
	int total;

	LinkedList();
	~LinkedList();
	void add(void*);
	void* get(int);
	void* remove(int);
	void** toArray(int& numberOfItems);
};
class LinkedListItem {
public:
	LinkedListItem* next;
	void* obj;
	LinkedListItem();
};

#endif
