# test_bag.py
# 
# Author: Seda den Boer
# Date: 01-02-2024
#
# Description: Test file for the Bag class.

import pytest
from classes.bag import Bag


def test_add():
    """
    Test the add() method of the Bag class by adding elements
    to a bag until it is full.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    assert test_bag.size == bag_capacity

    # Adding one more should raise an exception
    with pytest.raises(ValueError, match="Bag is full."):
        test_bag.add(bag_capacity + 1)


def test_remove():
    """
    Test the remove() method of the Bag class by adding elements
    to a bag and removing one of them.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    # Remove an element
    removed_index = 2
    test_bag.remove(removed_index)

    assert test_bag.size == bag_capacity - 1
    assert not test_bag.contains(removed_index)

    # Removing the same element again should raise an exception
    with pytest.raises(ValueError, match=f"Object with pool_index {removed_index} is not in the bag."):
        test_bag.remove(removed_index)


def test_pick():
    """
    Test the pick() method of the Bag class by adding elements
    to a bag and picking one of them.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Pick from an empty bag should raise an exception
    with pytest.raises(ValueError, match="Bag is empty."):
        test_bag.pick()

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    picked_index = test_bag.pick()

    assert picked_index is not None
    assert test_bag.contains(picked_index)


def test_contains():
    """
    Test the contains() method of the Bag class by adding elements
    to a bag and checking if they are in the bag.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    # Check if elements are in the bag
    for i in range(bag_capacity):
        assert test_bag.contains(i)

    # Check for an element not in the bag
    assert not test_bag.contains(bag_capacity + 1)


def test_get_number_occupied():
    """
    Test the get_number_occupied() method of the Bag class.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    assert test_bag.get_number_occupied() == 0

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    assert test_bag.get_number_occupied() == bag_capacity


def test_log(capsys):
    """
    Test the log() method of the Bag class by
    checking if the output is correct.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    test_bag.log()

    captured = capsys.readouterr()
    assert "elements" in captured.out
    assert "--" in captured.out


def test_str():
    """
    Test the __str__() method of the Bag class by
    checking if the output is correct.
    """
    bag_capacity = 5
    test_bag = Bag(pool_capacity=bag_capacity)

    # Add elements to the bag
    for i in range(bag_capacity):
        test_bag.add(i)

    assert str(test_bag) == str(test_bag.elements)