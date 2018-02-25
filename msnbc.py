import itertools
import numpy as np
import matplotlib.pyplot as plt
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file


class SeqToList:

    """Creates a SeqToList object from a .seq input file which can be converted into a list of lists"""

    def __init__(self, path):

        """path: str
               The path to the .seq file """

        self.path = path

    def seq_to_list_of_lists(self):

        """Takes a path to a .seq file and returns the source data as a list of lists"""

        f = open(self.path, 'r')
        lines = f.readlines()

        # remove ' \n' characters and parse the strings into a list of integers to be added to data
        data = []
        for string in lines:
            list_ints = string.replace(' \n', '').split(' ')
            ints = list(map(int, list_ints))
            data.append(ints)

        return data


class Markov:

    """Calculates and plots various quantities from a set of sequences of elements
    belonging to a finite set using Markov chains"""

    def __init__(self, data):

        """Creates a Markov object from list data"""

        self.data = data
        self.categories = set(itertools.chain(*[list(x) for x in data]))
        self.cardinality = len(self.categories)

        """data: list
               Source data containing sequences stored as a list of lists"""

    def get_transient_counts_matrix(self):

        """Returns a matrix containing the total counts of the number of transitions between
        elements in all sequences in the data. eg. a sequence with [..., 2, 3,....] > adds
         1 to the matrix element of column 2, row 3"""

        transient_counts_matrix = np.zeros(shape=(self.cardinality,
                                                  self.cardinality), dtype=int)
        for line in self.data:
            for i, j in zip(line[1:], line[:-1]):
                transient_counts_matrix[i-1][j-1] += 1

        return transient_counts_matrix

    def get_absorbing_counts_matrix(self):

        """Adds an additional element to the set representing the sequence 'termination' state
        i.e. the absorbing Markov state (the user going offline), and then appends this
        element to the end of each sequence. Then returns a matrix containing the total counts
        of the number of transitions between elements in all sequences in the data, including
        terminations"""

        absorbing_counts_matrix = np.zeros(shape=(self.cardinality + 1,
                                                  self.cardinality + 1), dtype=int)
        for line in data:
            line.append(0)
            for i, j in zip(line[1:], line[:-1]):
                absorbing_counts_matrix[i-1][j-1] += 1

        absorbing_counts_matrix[self.cardinality, self.cardinality] = 1

        return absorbing_counts_matrix

    @staticmethod
    def normalise(matrix):

        """matrix: An input matrix or array of arbitrary dimension.
        Returns the matrix calculated by dividing every element in a given column by the sum
        of all elements in that column"""

        column_sums = np.sum(matrix, axis=0)
        transition_matrix = matrix / column_sums[:]

        return transition_matrix

    def get_transient_markov_matrix(self):

        """Normalises by column the matrix returned by get_transient_counts_matrix. The columns
        add up to 1 and can be interpreted as probability mass functions"""

        m_t = self.normalise(self.get_transient_counts_matrix())

        return m_t

    def get_absorbing_markov_matrix(self):

        """Normalises by column the matrix returned by get_absorbing_counts_matrix. The columns
        add up to 1 and can be interpreted as probability mass functions"""

        m_a = self.normalise(self.get_absorbing_counts_matrix())

        return m_a

    def show_barchart(self, xmin=0, xmax=250, nbins=250, log='False'):

        """Produces a barchart showing the number of users (dependent variable) who visited a given
        number of pages (independent variable).

        xmin: int. Minimum no. pages visited
        xmax: int. Maximum no. pages visited
        nbins: int. Number of histogram bins
        log: bool. Use 'True' to view the dependent variable on a log10 scale"""

        # count the number of page visits for each sequence in data
        num_visits = []
        for i in range(len(self.data)):
            num_visits.append(len(self.data[i]))

        hist, edges = np.histogram(num_visits, bins=np.linspace(xmin, xmax, nbins))

        if log == 'True':

            p = figure(title="Log Barchart of Number of Page Views ({0} to {1} views)".format(xmin, xmax),
                       x_axis_label='Number of Page Views',
                       y_axis_label='Number of Users (log_10)',
                       tools="save"
                       )

            p.quad(top=np.log10(hist), bottom=0, left=edges[:-1], right=edges[1:], line_color="black")
            output_file('log_barchart.html', title="Log barchart")

        elif log == 'False':

            p = figure(title="Linear Barchart of Number of Page Views ({0} to {1} views)".format(xmin, xmax),
                       x_axis_label='Number of Page Views',
                       y_axis_label='Number of Users',
                       tools="save"
                       )

            p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="black")
            output_file('linear_barchart.html', title="Linear barchart")

        show(gridplot(p, ncols=1, plot_width=960, plot_height=540))

    def show_variability_scatterplot(self):

        """Produces a scatterplot with the proportion of the available page categories visited
        vs. the log10(total number of pages visited) in each sequence"""

        x = []
        y = []

        for line in self.data:
            y.append((len(np.unique(line)) / (len(self.categories))))
            x.append(np.log10(len(line)))

        plt.scatter(x, y, marker='.', s=1, alpha=0.4, c='black')
        plt.xlabel('Log_10 of number of page visits')
        plt.ylabel('Proportion of total page categories visited')
        plt.savefig('variability.png', dpi=384)
        plt.clf()

    def show_integer_ratios_scatterplot(self, page_category):

        """ Produces a scatterplot showing the ratio of unique page categories and total page views
        in each sequence vs. the proportion of total visits belonging to a particular page category

        page_category: int from 0 to 17 corresponding to the various categories
        """

        x = []
        y = []
        for line in self.data:
            x.append((len(np.unique(line)) / (len(line))))
            y.append(line.count(page_category) / len(line))

        plt.scatter(x, y, marker='.', s=1, alpha=0.25, c='black')
        plt.xlabel('Proportion of total page categories visited')
        plt.ylabel('Proportion of pages visited belonging to page '
                   'category {}'.format(page_category))
        plt.savefig('integer_ratios_{}.png'.format(page_category), dpi=384)
        plt.clf()

    def show_transition_matrix_heatmap(self):

        """Produces a heatmap image colour-coded to represent the probability of transition to a page
         category (row) given that the user is currently on a page category (column)"""

        fig, ax = plt.subplots()

        plt.imshow(self.get_absorbing_markov_matrix())
        plt.xticks(np.arange(0, len(self.categories) + 1, 1))
        plt.yticks(np.arange(0, len(self.categories) + 1, 1))

        category_labs = {0: "frontpage", 1: "news", 2: "tech", 3: "local",
                         4: "opinion", 5: "on-air", 6: "misc", 7: "weather",
                         8: "health", 9: "living", 10: "business", 11: "sports",
                         12: "summary", 13: "bbs", 14: "travel", 15: "msn-news",
                         16: "msn-sports", 17: "offline"}

        ax.set_xticklabels(list(category_labs.values()), rotation=90)
        ax.set_yticklabels(list(category_labs.values()))

        plt.gcf().subplots_adjust(bottom=0.2)
        plt.title('Probability that user will navigate to <choose row> \n given that user '
                  'is on <choose column>')
        plt.colorbar()
        plt.savefig('Transition probabilities matrix heatmap', dpi=384)
        plt.clf()

    def transition_probability_barchart(self, num_plots, input_matrix, category_num, savefigs='False'):

        """Displays bar charts showing the probability of the user navigating to each page category
        given that they started on a specified page category

        num_plots: Number of page-to-page transitions to consider (recommend < 20)
        input_matrix: Any column-normalised transition matrix, such as from get_absorbing_markov_matrix
        category_num: int from 0 to 17 corresponding to the various categories
        savefigs: Save .png files of probability barcharts. Default 'False'
        """

        # Apply transition matrix successively on initial state vector v

        identity = np.identity(len(self.categories) + 1)
        category_nums = np.arange(len(self.categories) + 1)

        category_labs = {0: "frontpage", 1: "news", 2: "tech", 3: "local",
                         4: "opinion", 5: "on-air", 6: "misc", 7: "weather",
                         8: "health", 9: "living", 10: "business", 11: "sports",
                         12: "summary", 13: "bbs", 14: "travel", 15: "msn-news",
                         16: "msn-sports", 17: "offline"}

        # initial vector corresponding to category of first page visited
        v = identity[category_num]

        # loop over successive transitions and update barchart. Save images optionally.
        plt.ion()
        for k in range(num_plots):
            plt.bar(category_nums, v, color='black', tick_label=category_labs.values())
            plt.ylim((0, 1))
            plt.title('Page navigation probability after {0} page \n'
                      'transitions given that first page is {1}'.format(k, category_labs[category_num]))
            plt.xticks(rotation='vertical')
            plt.ylabel('probability')
            plt.gcf().subplots_adjust(bottom=0.2)

            if savefigs == 'True':
                plt.savefig('{}_{}.png'.format(k, category_labs[category_num]))
            else:
                continue

            plt.pause(1)
            plt.clf()
            v = np.dot(input_matrix, v)
        plt.pause(2)
