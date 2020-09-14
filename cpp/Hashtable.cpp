
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


#include "Hashtable.h"

LinkedListItem::LinkedListItem() {
	next = NULL;
	obj = NULL;
}
LinkedList::LinkedList() {
	firstItem = NULL;
	lastItem = NULL;
	total = 0;
}
LinkedList::~LinkedList() {
	while (firstItem != NULL) {
		LinkedListItem* next = firstItem->next;
		delete firstItem;
		firstItem = next;
	}
	lastItem = NULL;
	total = 0;
}
void LinkedList::add(void* obj) {
	LinkedListItem* item = new LinkedListItem();
	item->obj = obj;
	if (lastItem == NULL) {
		firstItem = item;
	} else {
		lastItem->next = item;
	}
	lastItem = item;
	total++;
}
void* LinkedList::get(int index) {
	LinkedListItem* item = firstItem;
	for (int i=1;i<index;i++) {
		if (item == NULL) return NULL;
		item = item->next;
	}
	if (item == NULL) return NULL;
	return item->obj;
}
void* LinkedList::remove(int index) {
	if (firstItem == NULL) return NULL;

	LinkedListItem* item = firstItem;
	LinkedListItem* last = NULL;
	for (int i=1;i<index;i++) {
		if (item == NULL) return NULL;
		last = item;
		item = item->next;
	}
	if (item == NULL) return NULL;
	if (last == NULL) {//first item;
		firstItem = item->next;
	}
	if (item->next == NULL) {
		lastItem = NULL;
		if (last != NULL) {
			lastItem = last;
			last->next = NULL;
		}
	} else {
		if (last != NULL) {
			last->next = item->next;
		}
	}
	void *obj = item->obj;
	delete item;
	total--;
	return obj;
}
void** LinkedList::toArray(int &numberOfItems) {
	if (total == 0) return NULL;
	void** array = new void*[total];
	for (int i=0;i<total;i++) array[i] = NULL;
	numberOfItems = 0;
	LinkedListItem* item = firstItem;
	while (item != NULL) {
		array[numberOfItems++] = item->obj;
		item = item->next;
	}
	return array;
}


