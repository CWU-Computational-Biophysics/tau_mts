# plotlabelmanager.py

# Classes to help with labeling Matplotlib plots.

class PlotLabel:

    def __init__(self, key: object, text: str, plt_args: dict):
        # ensure that plt_args is at least an empty dict
        if plt_args is None:
            raise ValueError("plt_args cannot be None")
        if not isinstance(plt_args, dict):
            raise ValueError("plt_args must be a dict")

        # ensure text is a string
        if not isinstance(text, str):
            raise ValueError("text must be a string")

        # save data
        self.key = key
        self.text = text
        self.plt_args = plt_args

    def get_key(self) -> object:
        return self.key

    def get_text(self) -> str:
        return self.text

    def get_plt_args(self) -> dict:
        return self.plt_args

    def __repr__(self):
        return f"PlotLabel(key={self.key}, text='{self.text}', plt_args={self.plt_args})"



class PlotLabelManager:

    def __init__(self, args_map_dict: dict[str, list]):
        # define a dict to hold PlotLabel objects and access booleans
        self.labels = {}
        self.access = {}

        # store the argument map with lists reversed for popping
        self.arg_map = {k:v[::-1] for k, v in args_map_dict.items()}

        # check for inconsistent list lengths
        list_lens = set(len(v) for v in self.arg_map.values())
        if len(list_lens) > 1:
            # get the shortest list
            min_len = min(list_lens)

            # print warning
            print(f"Warning: not all lists in arg_map are the same length. Shortest list has {min_len} elements.")


    def get_plot_label(self, key: object) -> PlotLabel:
        # check if key exists
        if not self.key_exists(key):
            raise ValueError(f"Key '{key}' not found")

        return self.labels[key]


    def key_exists(self, key: object) -> bool:
        return key in self.labels


    def add(self, key: object, text: str, plt_args: dict = None) -> None:
        # check if key is already in list
        if self.key_exists(key):
            raise ValueError(f"Key '{key}' already exists")

        # if no plt_args are supplied, extract them from the arg_map
        if plt_args is None:
            # replace plt_args with an empty dict
            plt_args = {}

            # iterate over each key:list pair of the arg_map
            for arg_key, arg_list in self.arg_map.items():
                # check that arg_list is not empty
                if len(arg_list) == 0:
                    raise IndexError(f"arg_list for '{arg_key}' has been depleted of unique values")

                # the key of arg_key is the argument key in matplotlib.plt
                plt_args[arg_key] = arg_list.pop()

        # create a new PlotLabel object and append
        plot_label = PlotLabel(key, text, plt_args)
        self.labels[plot_label.get_key()] = plot_label
        self.access[plot_label.get_key()] = False


    def try_add(self, key: object, text: str, plt_args: dict = None) -> None:
        # check if key is not already in list
        if not self.key_exists(key):
            self.add(key=key, text=text, plt_args=plt_args)


    def get_args(self, key: object) -> dict:
        # prepare the return dict
        plot_label = self.get_plot_label(key)
        return_dict = plot_label.get_plt_args()
        return_dict['label'] = plot_label.get_text()

        # deduplicate the label entry for plotting
        if self.access[key]:
            return_dict['label'] = None
        else:
            self.access[key] = True

        return return_dict


    def __repr__(self):
        head = "PlotLabelManager:\n - "
        labels = "\n - ".join([f"{plot_label}" for plot_label in self.labels.values()])
        return_str = head + labels
        return return_str
