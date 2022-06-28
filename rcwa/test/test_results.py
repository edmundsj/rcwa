import pytest
from rcwa import Results
from matplotlib.figure import Axes, Figure
from matplotlib import pyplot as plt


@pytest.fixture()
def results():
    results_dict = {'rx': 1}
    res = Results(results_dict)
    yield res


def test_getitem(results):
    assert results['rx'] == 1


def test_setitem(results):
    with pytest.raises(Exception) as ex:
        results['rx'] = 2


def test_keys(results):
    assert results.keys() == {'rx': 1}.keys()


def test_values(results):
    assert list(results.values()) == list({'rx': 1}.values())


def test_items(results):
    assert results.items() == {'rx': 1}.items()


def test_plot_fig(results):
    fig_desired = Figure()
    fig_actual, ax_actual = results.plot(x='rx', y='rx', fig=fig_desired, ax=None)
    assert fig_actual is fig_desired
    assert isinstance(ax_actual, Axes)
    assert fig_actual.axes[0] is ax_actual


def test_plot_nofig(results):
    fig, ax = plt.subplots()
    fig_actual, ax_actual = results.plot(x='rx', y='rx', fig=None, ax=ax)
    assert isinstance(ax_actual, Axes)
    assert ax is ax_actual
    assert fig_actual is None
    assert fig.axes[0] is ax_actual


def test_plot_nofig_noax(results):
    fig_actual, ax_actual = results.plot(x='rx', y='rx', fig=None, ax=None)
    assert isinstance(fig_actual, Figure)
    assert isinstance(ax_actual, Axes)
    assert fig_actual.axes[0] is ax_actual



