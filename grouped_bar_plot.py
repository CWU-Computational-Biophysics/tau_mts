# grouped_bar_plot.py
# functions used for making categorical grouped bar plot

# imports
import numpy as np
from matplotlib import colormaps


# define a function to make a grouped bar plot
def grouped_bar_plot(fig, ax,
                     group_data: dict,
                     group_labels: list,
                     colors: list = None,
                     bar_width: float = 0.15,
                     multiplier: float = 0,
                     **kwargs):
    """Plots a grouped bar plot

    Arguments:
        fig -- _description_
        ax -- _description_
        group_data -- _description_

    Keyword Arguments:
        colors -- _description_ (default: {None})
        bar_width -- _description_ (default: {0.2})
        multiplier -- _description_ (default: {0})

    Raises:
        ValueError: _description_

    Returns:
        _description_
    """

    # extract groups from keys of group_data
    group_count = len(group_labels)

    # the default value of colors should be a built in color list
    if colors is None:
        colors = colormaps["tab10"].colors
    else:
        # check if the length of colors is greater than the number of groups
        if len(colors) < group_count:
            raise ValueError("The number of colors provided is less than the number of groups")

    # define the label locations
    x_locations = np.arange(group_count)

    # plot the bars
    for index, (attribute, measurement) in enumerate(group_data.items()):
        # calculate the offset location of each bar
        offset = bar_width * multiplier

        # draw the bar and insert kwargs
        ax.bar(
            x_locations + offset,
            measurement,
            bar_width,
            label=attribute,
            color=colors[index],
            **kwargs
        )

        # iterate the multiplier to space the next bar
        multiplier += 1

    # set the x ticks
    ax.set_xticks(x_locations + bar_width, group_labels)

    return fig, ax
