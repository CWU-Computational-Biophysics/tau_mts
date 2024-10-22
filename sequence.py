# sequence.py

# A class for sequence information.
# TODO: Add docstrings


class Sequence:

    def __init__(self, start_index: int, end_index: int, type: int, grid_dict: dict):
        # check that protein type is in the values of grid_dict
        if type not in grid_dict.values():
            raise ValueError(f"Invalid protein type: {type}")

        # save the indices
        self.start_index = start_index
        self.end_index = end_index
        self.length = self.end_index - self.start_index

        # ensure length is positive
        if self.length <= 0:
            raise ValueError(f"Invalid start or end index: {start_index}, {end_index}")

        # save the protein type and name
        self.type = type
        self.protein = list(grid_dict.keys())[list(grid_dict.values()).index(self.type)]


    def __str__(self):
        return f"Sequence of {self.length} {self.protein} starting at {self.start_index}"


    def __repr__(self):
        return self.__str__()


    def get_start_index(self) -> int:
        return self.start_index


    def get_end_index(self) -> int:
        return self.end_index


    def get_length(self) -> int:
        return self.length


    def get_type(self) -> int:
        return self.type
