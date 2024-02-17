# test_pool.py
# 
# Author: Seda den Boer
# Date: 01-02-2024
#
# Description: Test file for the Pool class.

import pytest
from classes.pool import Pool
from classes.triangle import Triangle


def test_occupy():
    """
    Test the occupy() function of the Pool class by
    occupying all spaces in the pool.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    assert test_pool.get_number_occupied() == pool_capacity

    # Occupying one more should raise an exception
    with pytest.raises(Exception, match="Pool is full."):
        test_pool.occupy(Triangle())


def test_free():
    """
    Test the free() function of the Pool class by
    freeing all spaces in the pool.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    # Free an index
    freed_index = occupied_indices[2]
    test_pool.free(freed_index)

    assert test_pool.get_number_occupied() == pool_capacity - 1
    assert not test_pool.contains(freed_index)

    # Freeing the same index again should raise an exception
    with pytest.raises(Exception, match="Index is not valid."):
        test_pool.free(freed_index)


def test_get():
    """
    Test the get() function of the Pool class by getting all
    objects from the pool.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    # Get objects from the pool
    for index in occupied_indices:
        assert isinstance(test_pool.get(index), Triangle)


def test_contains():
    """
    Test the contains() function of the Pool class by
    checking if all occupied indices are in the pool.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    # Check if indices are in the pool
    for index in range(pool_capacity):
        if index in occupied_indices:
            assert test_pool.contains(index)
        else:
            assert not test_pool.contains(index)

    # Check for an index not in the pool
    assert not test_pool.contains(pool_capacity + 1)


def test_get_number_occupied():
    """
    Test the get_number_occupied() function of the Pool class.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    assert test_pool.get_number_occupied() == 0

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    assert test_pool.get_number_occupied() == pool_capacity


def test_log(capsys):
    """
    Test the log() function of the Pool class by
    checking if the output is correct.
    """
    pool_capacity = 5
    test_pool = Pool(capacity=pool_capacity)

    # Occupy spaces in the pool
    occupied_indices = [test_pool.occupy(Triangle()) for _ in range(pool_capacity)]

    test_pool.log()

    captured = capsys.readouterr()
    assert "elements" in captured.out
    assert "x" in captured.out
    assert "p" in captured.out
    assert "used indices" in captured.out
    assert "first" in captured.out