import matplotlib.pyplot as plt


class Plot:
    def __init__(self):
        self.fig = None
        self.ax = None
        self.labels = list()

    def plot(self, x, y, label: str = None, num="simulated_PDF"):
        if self.fig:
            if num !=self.fig._label:
                self.fig, self.ax = plt.subplots(1, 1, num=num)
        else:
            self.fig, self.ax = plt.subplots(1, 1, num=num)

        if label:
            self.labels.append(label)
            self.ax.plot(x, y, label=label)
        else:
            self.ax.plot(x, y)

        self.finalize()


    def finalize(self):
        self.ax.set_xlabel(r"$r (\AA)$")
        self.ax.set_ylabel(r"$G (\AA^{-2})$")
        # if self.labels:
        #     self.ax.legend()

    def show(self):
        if self.fig:
            self.fig.canvas.draw_idle()
            self.fig.canvas.flush_events()
        else:
            plt.show()
