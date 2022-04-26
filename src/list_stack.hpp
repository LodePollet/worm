
#pragma once
#include <algorithm>
#include <list>
#include <vector>
#include <utility>

template <typename T>
class list_stack {
public:
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;
    typedef const T& const_reference;
    typedef const T* const_pointer;
    typedef std::pair<T, size_t> pair_type;
    typedef std::list<pair_type> list_type;
    typedef typename list_type::size_type size_type;
    class iterator : public list_type::iterator {
    public:
        friend list_stack;
        typedef T value_type;
        typedef T& reference;
        typedef T* pointer;
        typedef std::bidirectional_iterator_tag iterator_category;
        iterator() : list_type::iterator() {}
        iterator(const typename list_type::iterator& it)
            : list_type::iterator(it) {}
        inline reference operator*() const {
            return list_type::iterator::operator*().first;
        }
        inline pointer operator->() const {
            return &(list_type::iterator::operator*().first);
        }
    private:
        inline size_t& stack_pos() {
            return list_type::iterator::operator*().second;
        }
    };
    class const_iterator : public list_type::const_iterator {
    public:
        friend list_stack;
        typedef T value_type;
        typedef const T& const_reference;
        typedef const T* const_pointer;
        typedef std::bidirectional_iterator_tag iterator_category;
        const_iterator() : list_type::const_iterator() {}
        const_iterator(const typename list_type::iterator& it)
            : list_type::const_iterator(it) {}
        const_iterator(const typename list_type::const_iterator& it)
            : list_type::const_iterator(it) {}
        inline const_reference operator*() const {
            return list_type::const_iterator::operator*().first;
        }
        inline const_pointer operator->() const {
            return &(list_type::const_iterator::operator*().first);
        }
    private:
        inline const size_t& stack_pos() const {
            return list_type::const_iterator::operator*().second;
        }
    };

    typedef std::vector<iterator> stack_type;

    iterator insert(const_iterator pos, const T& val) {
        iterator it = list.insert(pos, std::make_pair(val, stack.size()));
        stack.push_back(it);
        return it;
    }

    iterator insert(const_iterator pos, T&& val) {
        iterator it = list.insert(pos,
                                  std::make_pair(std::move(val), stack.size()));
        stack.push_back(it);
        return it;
    }

    iterator erase(const_iterator pos) {
        iterator last = stack.back();
        stack.pop_back();
        last.stack_pos() = pos.stack_pos();
        stack[pos.stack_pos()] = last;
        return list.erase(pos);
    }

    iterator at(size_t i) {
        return stack[i];
    }

    const_iterator at(size_t i) const {
        return stack[i];
    }

    iterator begin() {
        return list.begin();
    }

    const_iterator begin() const {
        return list.begin();
    }

    iterator end() {
        return list.end();
    }

    const_iterator end() const {
        return list.end();
    }

    reference front() {
        return list.front().first;
    }

    const_reference front() const {
        return list.front().first;
    }

    reference back() {
        return list.back().first;
    }

    const_reference back() const {
        return list.back().first;
    }

    size_type size() const {
        return list.size();
    }

    void clear() {
        stack.clear();
        list.clear();
    }

    void sort() {
        list.sort();
    }

    template <class Compare>
    void sort(Compare cmp) {
        list.sort(cmp);
    }

private:
    list_type list;
    stack_type stack;
};
