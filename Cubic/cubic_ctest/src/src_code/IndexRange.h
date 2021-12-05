/*
 * This piece of code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <assert.h>

class IndexRange {
private:
	int64_t start_ = 0;
	int64_t size_ = 0;

public:
	constexpr IndexRange() = default;

	constexpr explicit IndexRange(int64_t size) : start_(0), size_(size)
	{
		assert(size >= 0);
	}

	constexpr IndexRange(int64_t start, int64_t size) : start_(start), size_(size)
	{
		assert(start >= 0);
		assert(size >= 0);
	}

	class Iterator {
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = int64_t;
		using pointer = const int64_t*;
		using reference = const int64_t&;
		using difference_type = std::ptrdiff_t;

	private:
		int64_t current_;

	public:
		constexpr explicit Iterator(int64_t current) : current_(current)
		{
		}

		constexpr Iterator& operator++()
		{
			current_++;
			return *this;
		}

		constexpr Iterator operator++(int) const
		{
			Iterator copied_iterator = *this;
			++copied_iterator;
			return copied_iterator;
		}

		constexpr friend bool operator!=(const Iterator& a, const Iterator& b)
		{
			return a.current_ != b.current_;
		}

		constexpr friend bool operator==(const Iterator& a, const Iterator& b)
		{
			return a.current_ == b.current_;
		}

		constexpr int64_t operator*() const
		{
			return current_;
		}
	};

	constexpr Iterator begin() const
	{
		return Iterator(start_);
	}

	constexpr Iterator end() const
	{
		return Iterator(start_ + size_);
	}

	/**
	 * Access an element in the range.
	 */
	constexpr int64_t operator[](int64_t index) const
	{
		assert(index >= 0);
		assert(index < this->size());
		return start_ + index;
	}

	/**
	 * Two ranges compare equal when they contain the same numbers.
	 */
	constexpr friend bool operator==(IndexRange a, IndexRange b)
	{
		return (a.size_ == b.size_) && (a.start_ == b.start_ || a.size_ == 0);
	}

	/**
	 * Get the amount of numbers in the range.
	 */
	constexpr int64_t size() const
	{
		return size_;
	}

	/**
	 * Create a new range starting at the end of the current one.
	 */
	constexpr IndexRange after(int64_t n) const
	{
		assert(n >= 0);
		return IndexRange(start_ + size_, n);
	}

	/**
	 * Create a new range that ends at the start of the current one.
	 */
	constexpr IndexRange before(int64_t n) const
	{
		assert(n >= 0);
		return IndexRange(start_ - n, n);
	}

	/**
	 * Get the first element in the range.
	 * Asserts when the range is empty.
	 */
	constexpr int64_t first() const
	{
		assert(this->size() > 0);
		return start_;
	}

	/**
	 * Get the last element in the range.
	 * Asserts when the range is empty.
	 */
	constexpr int64_t last() const
	{
		assert(this->size() > 0);
		return start_ + size_ - 1;
	}

	/**
	 * Get the element one after the end. The returned value is undefined when the range is empty.
	 */
	constexpr int64_t one_after_last() const
	{
		return start_ + size_;
	}

	/**
	 * Get the first element in the range. The returned value is undefined when the range is empty.
	 */
	constexpr int64_t start() const
	{
		return start_;
	}

	/**
	 * Returns true when the range contains a certain number, otherwise false.
	 */
	constexpr bool contains(int64_t value) const
	{
		return value >= start_ && value < start_ + size_;
	}

	/**
	 * Returns a new range, that contains a sub-interval of the current one.
	 */
	constexpr IndexRange slice(int64_t start, int64_t size) const
	{
		assert(start >= 0);
		assert(size >= 0);
		int64_t new_start = start_ + start;
		assert(new_start + size <= start_ + size_ || size == 0);
		return IndexRange(new_start, size);
	}
	constexpr IndexRange slice(IndexRange range) const
	{
		return this->slice(range.start(), range.size());
	}


	friend std::ostream& operator<<(std::ostream& stream, IndexRange range)
	{
		stream << "[" << range.start() << ", " << range.one_after_last() << ")";
		return stream;
	}
};